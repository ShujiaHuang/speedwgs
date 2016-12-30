#!/usr/bin/env perl

use strict;
use warnings;
use POSIX;
use POSIX ":sys_wait_h";
use File::Basename;

#args
my $useDisk = shift @ARGV;

if (eof(STDIN)) {
    print STDERR "no input, exit now.\n";
    exit;
}

my $bin = "gene_bin.tar.gz";

use MIME::Base64;
my $delim = "\x01";

sub system2 {
    my $cmd = shift;
    my $t = time;
    print STDERR strftime("%F %T", localtime), " >> start command: $cmd\n";
    0 == system($cmd) or die "command failed: $cmd\n";
    $t = time - $t;
    print STDERR strftime("%F %T", localtime), " << command finished in $t seconds: $cmd\n";
}

my %chr_length = (
    "1" => 249250621,
    "2" => 243199373,
    "3" => 198022430,
    "4" => 191154276,
    "5" => 180915260,
    "6" => 171115067,
    "7" => 159138663,
    "8" => 146364022,
    "9" => 141213431,
    "10" => 135534747,
    "11" => 135006516,
    "12" => 133851895,
    "13" => 115169878,
    "14" => 107349540,
    "15" => 102531392,
    "16" => 90354753,
    "17" => 81195210,
    "18" => 78077248,
    "19" => 59128983,
    "20" => 63025520,
    "21" => 48129895,
    "22" => 51304566,
    "X" => 155270560,
    "Y" => 59373566,
    "MT" => 16569
    );

my $pid = fork();

if ($pid == 0) { # child process do the computing

    my $reducer_id = -1;
    my $tmp_root = -d "/var/run" ? "/var/run" : ".";
    if (defined($useDisk) and $useDisk eq "disk") {
        $tmp_root = ".";
    }
    my $code = 0; # subprocess exit code

    # format input and samtools rmdup, here bam will be written to disk
    my $bam = "$tmp_root/$$-rmdup.bam";
    open RMDUP, "| $bin/samtools rmdup - $bam 2>/dev/null";
    my $header_printed = 0;
    my @key_seq = ();
    my @seq2rname = (1 .. 22, "X", "Y", "MT");
    my %chr_seq;
    for (my $i = 0; $i <= $#seq2rname; $i++) {
        $chr_seq{$seq2rname[$i]} = $i + 1;
    }
    my %chr2interval;
    print STDERR strftime("%F %T", localtime), " >> start command: samtools rmdup\n";
    my $rmdup_time = time;
    while (<>) {
        my @fields = split /$delim/;
        if ($fields[1] == 0 and $fields[2] == 0) { # header
            if (!$header_printed) {
                $reducer_id = $fields[0];
                print STDERR "reducer $reducer_id\n";
                $header_printed = 1;
                print RMDUP decode_base64($fields[3]);
            }
        } else {
            print RMDUP $fields[3];
            my $rname = $seq2rname[$fields[1] - 1];
            if (exists($chr2interval{$rname})) {
                $chr2interval{$rname}[1] = $fields[2];
            } else {
                @{$chr2interval{$rname}} = ($fields[2], $fields[2]);
                @key_seq = (@key_seq, $rname);
            }
        }
    }
    close RMDUP;
    $code = $? >> 8;
    exit($code) unless ($code == 0);
    $rmdup_time = time - $rmdup_time;
    print STDERR strftime("%F %T", localtime), " << command finished in $rmdup_time seconds: samtools rmdup\n";

    sub get_low_bound {
        my $val = shift;
        my $low = $val - 500;
        if ($low < 1) {
            return(1);
        }
        return($low);
    }

    sub get_up_bound {
        my $chr = shift;
        my $val = shift;
        my $up = $val + 600;
        if ($up > $chr_length{$chr}) {
            return($chr_length{$chr});
        }
        return($up);
    }

    print STDERR strftime("%F %T", localtime), " >> calculate input interval\n";
    my $interval = "$tmp_root/$$-interval.list";
    open INTERVAL, ">", $interval;
    my $item;
    foreach $item (@key_seq) {
        my $str = $item.':'.get_low_bound($chr2interval{$item}[0]).'-'.get_up_bound($item, $chr2interval{$item}[1])."\n";
        print INTERVAL $str;
        print STDERR $str; # also print to stderr
    }
    close INTERVAL;
    print STDERR strftime("%F %T", localtime), " << calculate input interval\n";
    die "reducer id not determined.\n" if ($reducer_id < 0);

    # samtools index
    system2("$bin/samtools index $bam");

    # gatk RealignerTargetCreator
    my $bed = "$tmp_root/$$-target.bed";
    system2("[ -e alijre-1.6.0_45 ] || ln -s /usr/ali/jre alijre-1.6.0_45");
    system2("$bin/gatk --fix_misencoded_quality_scores -T RealignerTargetCreator -l OFF -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -L $interval -I $bam -o $bed");

    # gatk IndelRealigner and output
    # -l off : http://gatkforums.broadinstitute.org/wdl/discussion/3830/piping-gatk-output-to-stdout
    my $realignedBam = "$tmp_root/$$-realignedBam.bam";
    system2("$bin/gatk -T IndelRealigner -rf NotPrimaryAlignment -L $interval -l OFF -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -I $bam -known 1000G_phase1.indels.tar.gz/1000G_phase1.indels.b37.vcf -known Mills_and_1000G_gold_standard.indels.tar.gz/Mills_and_1000G_gold_standard.indels.b37.vcf -targetIntervals $bed -o $realignedBam");

    # gatk BaseRecalibrator
    my $baseRecReport = "$tmp_root/$$-base-recalibration-report.grp";
    system2("$bin/gatk -T BaseRecalibrator -L $interval -l OFF --disable_indel_quals -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -knownSites 1000G_phase1.indels.tar.gz/1000G_phase1.indels.b37.vcf -knownSites Mills_and_1000G_gold_standard.indels.tar.gz/Mills_and_1000G_gold_standard.indels.b37.vcf -knownSites dbsnp_138.tar.gz/dbsnp_138.b37.vcf -I $realignedBam -o $baseRecReport");

    # gatk PrintReads
    my $outBam = "$tmp_root/$$-recOutBam.bam";
    system2("$bin/gatk -T PrintReads -L $interval -l OFF -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -I $realignedBam -BQSR $baseRecReport -o $outBam");

    # gatk UnifiedGenotyper
    my $unifiedCmd = "$bin/gatk -T UnifiedGenotyper -L $interval -l OFF -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -I $outBam --dbsnp dbsnp_138.tar.gz/dbsnp_138.b37.vcf |";
    my $unifiedTime = time;
    print STDERR strftime("%F %T", localtime), " >> start command: $unifiedCmd\n";
    my $header = "";
    my $out_header_printed = 0;
    open OUT, $unifiedCmd;
    while (<OUT>) {
        if (m/^#/) {
            $header .= $_;
            next;
        }
        if (!$out_header_printed) {
            # set chr seq and pos to -1 for sort
            print join($delim, -1, -1, encode_base64($header, "")), "\n";
            $out_header_printed = 1;
        }
        my ($rname, $pos, $others) = split("\t", $_, 3);
        print join($delim, $chr_seq{$rname}, $pos, $_);
    }
    close OUT;
    $code = $? >> 8;
    exit($code) unless ($code == 0);
    $unifiedTime = time - $unifiedTime;
    print STDERR strftime("%F %T", localtime), " << command finished in $unifiedTime seconds: $unifiedCmd\n";

    # TODO: clean all tmp files?

    exit(0);

} elsif ($pid) { # parent process monitors and heart beats
    my $c = 0;
    while (1) {
        my $child = waitpid($pid, WNOHANG);
        if ($child != 0) {
            $? == 0 ? print STDERR "everything seems good, exit\n" : print STDERR "something bad happened.\n";
            exit($? >> 8);
        }
        if ($c++ > 60 * 5) { # default heartbreak timeout is 10 min
            print STDERR "heartbeat\n";
            $c = 0;
        }
        sleep(1);
    }
}
