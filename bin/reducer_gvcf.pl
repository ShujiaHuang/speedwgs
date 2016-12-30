#!/usr/bin/env perl
use strict;
use warnings;
use POSIX;
use POSIX ":sys_wait_h";
use File::Basename;
use Gene_util;
use MIME::Base64;


Gene_util::check_eof();
#args
my $tmp_storage = shift @ARGV;
my $gatk_script = shift @ARGV;
my $output_format = shift @ARGV;
my $input_table_name = shift @ARGV;

die "Usage: $0 tmp_storage gatk_script output_format input_table_name\n"
    unless (defined($tmp_storage) and defined($gatk_script) and defined($output_format) and defined($input_table_name));

my $bin = Gene_util::get_gene_bin();
my $delim = Gene_util::get_field_delim();
my %std_ref = Gene_util::get_std_ref();
my $gatk_cmd = "$bin/$gatk_script";

sub get_realigner_cmd {
    my $interval = shift;
    my $bam = shift;
    my $bed = shift;

    my $opt = $gatk_script eq "gatk3.5" ? "" : "--fix_misencoded_quality_scores";
    return "$gatk_cmd " . $opt . " -T RealignerTargetCreator -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -L $interval -I $bam -o $bed";
}

sub get_output_cmd {
    my $interval = shift;
    my $outBam = shift;
    my $out_file = shift;

    if ($output_format eq "vcf") {
        return "$gatk_cmd -T UnifiedGenotyper -L $interval -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -I $outBam --dbsnp dbsnp_138.tar.gz/dbsnp_138.b37.vcf -o $out_file";
    } elsif ($gatk_script eq "gatk3.5") { # $output_format = gvcf
        return "$gatk_cmd -T HaplotypeCaller -L $interval -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -I $outBam --emitRefConfidence GVCF -o $out_file";
    } else {
        die "Cannot execute with '$gatk_script' and '$output_format'\n";
    }
}

my $pid = fork();
if ($pid == 0) { # child process do the computing

    my $reducer_id = -1;
    my $tmp_root = -d "/var/run" ? "/var/run" : ".";
    if (defined($tmp_storage) and $tmp_storage eq "disk") {
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
    my $rmdup_time = time;
    while (<>) {
        my ($part, $chr_seq, $pos, $content) = split /$delim/;
        if ($chr_seq == 0 and $pos == 0) { # header
            if (!$header_printed) {
                $reducer_id = $part;
                print STDERR strftime("%F %T", localtime), " >> start worker: $reducer_id\n";
                print STDERR strftime("%F %T", localtime), " >> start command: samtools rmdup \n";
                $header_printed = 1;
                print RMDUP decode_base64($content);
            }
        } else {
            print RMDUP $content;
            my $rname = $seq2rname[$chr_seq - 1]; # record start end interval
            if (exists($chr2interval{$rname})) {
                $chr2interval{$rname}[1] = $pos;
            } else {
                @{$chr2interval{$rname}} = ($pos, $pos);
                @key_seq = (@key_seq, $rname);
            }
        }
    }
    close RMDUP;
    $code = $? >> 8;
    exit($code) unless ($code == 0);

    $rmdup_time = time - $rmdup_time;
    print STDERR strftime("%F %T", localtime), " << command finished in $rmdup_time seconds: rmdup\n";

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
        if ($up > $std_ref{$chr}) {
            return($std_ref{$chr});
        }
        return($up);
    }

    print STDERR strftime("%F %T", localtime), " >> calculate input interval\n";
    my $interval = "$tmp_root/$$-interval.list";
    open INTERVAL, ">", $interval;
    my $item;
    print STDERR "Real interval:\n";
    foreach $item (@key_seq) {
        my $s = $item.':'.get_low_bound($chr2interval{$item}[0]).'-'.get_up_bound($item, $chr2interval{$item}[1])."\n";
        print INTERVAL $s;
        print STDERR $s;
    }
    close INTERVAL;
    print STDERR strftime("%F %T", localtime), " << calculate input interval\n";
    die "reducer id not determined.\n" if ($reducer_id < 0);

    # samtools index
    Gene_util::system2("$bin/samtools index $bam");

    # gatk RealignerTargetCreator
    my $bed = "$tmp_root/$$-target.bed";
    Gene_util::system2("[ -e alijre-1.6.0_45 ] || ln -s /usr/ali/jre1.7.0_79 alijre-1.6.0_45");
    Gene_util::system2(&get_realigner_cmd($interval, $bam, $bed));

    # gatk IndelRealigner and output
    # -l off : http://gatkforums.broadinstitute.org/wdl/discussion/3830/piping-gatk-output-to-stdout
    my $realignedBam = "$tmp_root/$$-realignedBam.bam";
    Gene_util::system2("$gatk_cmd -T IndelRealigner -rf NotPrimaryAlignment -L $interval -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -I $bam -known 1000G_phase1.indels.tar.gz/1000G_phase1.indels.b37.vcf -known Mills_and_1000G_gold_standard.indels.tar.gz/Mills_and_1000G_gold_standard.indels.b37.vcf -targetIntervals $bed -o $realignedBam");
    unlink($bam, "$bam.bai", $bed);

    # gatk BaseRecalibrator
    my $baseRecReport = "$tmp_root/$$-base-recalibration-report.grp";
    Gene_util::system2("$gatk_cmd -T BaseRecalibrator -L $interval --disable_indel_quals -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -knownSites 1000G_phase1.indels.tar.gz/1000G_phase1.indels.b37.vcf -knownSites Mills_and_1000G_gold_standard.indels.tar.gz/Mills_and_1000G_gold_standard.indels.b37.vcf -knownSites dbsnp_138.tar.gz/dbsnp_138.b37.vcf -I $realignedBam -o $baseRecReport");

    # gatk PrintReads
    my $outBam = "$tmp_root/$$-recOutBam.bam";
    Gene_util::system2("$gatk_cmd -T PrintReads -L $interval -R human_g1k_v37.tar.gz/human_g1k_v37.fasta -I $realignedBam -BQSR $baseRecReport -o $outBam");
    unlink($realignedBam, substr($realignedBam, 0, -1)."i");
    unlink($baseRecReport);

    # gatk HaplotypeCaller, ouptut gvcf
    # XXX: have to make a named pipe and fork twice because this gatk3.5 bug:
    # http://gatkforums.broadinstitute.org/gatk/discussion/6950/haplotypecaller-runtime-error-with-emitrefconfidence-gvcf
    use POSIX qw(mkfifo);
    my $outGVCF = "$tmp_root/$$-$reducer_id.g.vcf";
    $code = mkfifo($outGVCF, 0700);
    die "FATAL: failed to create named pipe\n" unless defined($code);
    my $unifiedCmd = &get_output_cmd($interval, $outBam, $outGVCF);
    my $unifiedTime = time;
    print STDERR strftime("%F %T", localtime), " >> start command: $unifiedCmd\n";

    my $gatk_pid = fork();
    die "FATAL: failed to fork process for GATK HaplotypeCaller\n" unless defined($gatk_pid);
    if ($gatk_pid == 0) {
        exec $unifiedCmd;
    }
    my $pipe_pid = fork();
    die "FATAL: failed to fork process for pipe output\n" unless defined($pipe_pid);
    if ($pipe_pid == 0) {
        my $header = "";
        my $out_header_printed = 0;
        open GVCF, "<", $outGVCF or die "FATAL: can not open file $outGVCF\n";
        while (<GVCF>) {
            if (m/^#/) {
                $header .= $_;
                next;
            }
            if (!$out_header_printed) {
                # set chr seq and pos to -1 for sort
                print join($delim, $reducer_id, -1, -1, encode_base64($header, "")), "\n";
                $out_header_printed = 1;
            }
            my ($rname, $pos, $others) = split("\t", $_, 3);
            print join($delim, $reducer_id, $chr_seq{$rname}, $pos, $_);
        }
        close GVCF;
        exit 0;
    }
    waitpid($gatk_pid, 0);
    $code = $? >> 8;
    exit($code) unless ($code == 0);
    waitpid($pipe_pid, 0);
    $code = $? >> 8;
    exit($code) unless ($code == 0);
    $unifiedTime = time - $unifiedTime;
    print STDERR strftime("%F %T", localtime), " << command finished in $unifiedTime seconds: $unifiedCmd\n";
    unlink $outGVCF, $outBam, substr($outBam, 0, -1)."i", $interval;
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
