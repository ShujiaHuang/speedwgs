#!/usr/bin/env perl

use strict;
use warnings;
use POSIX;
use POSIX ":sys_wait_h";
use File::Basename;

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

my $pid = fork();

if ($pid == 0) { # child process do the computing

    my $reducer_id = -1;

    # format input and samtools rmdup, here bam will be written to disk
    my $bam = "tmp.bam";
    open REDUP, "| $bin/samtools rmdup - $bam";
    my $header_printed = 0;
    while (<>) {
        my @fields = split /$delim/;
        if ($fields[1] == 0 and $fields[2] == 0) { # header
            if (!$header_printed) {
                $reducer_id = $fields[0];
                $header_printed = 1;
                print REDUP decode_base64($fields[3]);
            }
        } else {
            print REDUP $fields[3];
        }
    }
    close REDUP;
    die "reducer id not determined.\n" if ($reducer_id < 0);

    # samtools index
    system2("$bin/samtools index $bam");

    # gatk RealignerTargetCreator
    my $bed = "tmp.bed";
    my $interval = "interval.list";
    system2("$bin/samtools view $bam | $bin/get_sam_interval.py > $interval");

    open IN, "<", $interval;
    my $outputStr = $reducer_id.$delim;
    my $line;
    while (!eof(IN)) {
        chomp($line=<IN>);
        $outputStr = $outputStr.$line.",";
    }
    print $outputStr."\n";
    exit(0);

} elsif ($pid) { # parent process monitors and heart beats
    my $c = 0;
    while (1) {
        my $child = waitpid($pid, WNOHANG);
        if ($child != 0) {
            $? == 0 ? print STDERR "everything seems good, exit\n" : print STDERR "something bad happened.\n";
            exit($?);
        }
        if ($c++ > 60 * 5) { # default heartbreak timeout is 10 min
            print STDERR "heartbeat\n";
            $c = 0;
        }
        sleep(1);
    }
}
