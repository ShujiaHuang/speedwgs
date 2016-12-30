#!/usr/bin/env perl

use strict;
use warnings;
use POSIX;
use File::Basename;
use threads;
use Thread::Queue;
use Getopt::Long;

my $conf_file;
my $thread_num = 4;
my $overwrite = 0;

sub print_usage {
    print <<EOF;
usage: $0 [options] list_file ...
options:
    -c --config     ODPS config file
    -t --thread     number of threads to use, default $thread_num
    -o --overwrite  overwrite data in existing table, default disabled
    -h --help       print this usage
EOF
;
}

GetOptions("c|config=s" => \$conf_file,
           "t|thread=i" => \$thread_num,
           "o|overwrite" => \$overwrite,
           "h|help" => sub { &print_usage(); exit 0; })
    or exit 1;

die "config file not found\n" unless (defined($conf_file) and -e $conf_file);
die "thread number must >= 1: $thread_num\n" unless ($thread_num >= 1);

# make odpscmd/tattoo in PATH
my $bin_dir = dirname($0);
my $path= $ENV{'PATH'};
$ENV{'PATH'} = "$path:$bin_dir/odpscmd/bin:$bin_dir/tattoo/bin:$bin_dir/data_upload";
# lead all child processes
setsid();

sub debug {
    my $msg = shift;
    print STDERR strftime("%F %T", localtime), " $msg\n";
}

my @failures;
sub fail {
    my $msg = shift;
    push @failures, $msg;
    debug("XX $msg");
}

# return 0(false) or 1(true)
# TODO: calc md5 might be slow if hard disk is broken, is timeout monitor necessary?
sub check_md5 {
    my $file = shift;
    my $md5 = shift;

    my $ret = `md5sum $file | head -c 32`;
    if ($? != 0) {
        fail("failed to calculate md5sum for $file");
        return 0;
    }
    if ($ret ne $md5) {
        fail("md5sum not match: $file: $ret vs $md5");
        return 0;
    }
    return 1;
}

my $delim = "\x00";
my $q = Thread::Queue->new();
my @runners;
for (1..$thread_num) {
    push @runners, async {
        while (my $item = $q->dequeue()) {
            my ($table, $fq1, $fq2, $pos, $adpt1, $adpt2) = split /$delim/, $item;
            # zhuyun's upload terminal has 100~200MBytes/s bandwidth, so tattoo use 8 threads to upload
            my $cmd = "merge_fastq.pl $fq1 $fq2 $table $conf_file 8 '$adpt1' '$adpt2'";
            debug(">> $pos - upload to $table: $fq1 + $fq2");
            if (0 == system($cmd)) {
                debug("<< $pos - upload complete");
            } else {
                fail("$pos - upload failed: $fq1 + $fq2");
            }
        }
    };
}

my %tables;
sub upload {
    my $table = shift;
    my $fq1 = shift;
    my $fq2 = shift;
    my $pos = shift;
    my $adpt1 = shift;
    my $adpt2 = shift;
    unless (exists($tables{$table})) {
        my $cmd_truncate = $overwrite ? "truncate table $table" : "";
        my $cmd = "odpscmd --config=$conf_file -e 'create table if not exists $table (id1 string, seq1 string, opt_id1 string, quality1 string, id2 string, seq2 string, opt_id2 string, quality2 string);$cmd_truncate'";
        debug($overwrite ? "-- create and truncate table $table" : "-- create table $table");
        unless (0 == system($cmd)) {
            fail("$pos: failed to create table $table");
            return;
        }
        $tables{$table}++;
    }
    $q->enqueue(join($delim, $table, $fq1, $fq2, $pos, $adpt1, $adpt2));
}

for my $l (@ARGV) {
    debug(">> $l: read list file");
    my $f;
    unless (open $f, "<", $l) {
        fail("$l: can not open file");
        next;
    }
    my $p = dirname($l);
    while (<$f>) {
        next if m/^\s*$/;  # ignore empty line
        next if m/^\s*#/;  # ignore comment line
        s/[\r\n]+$//;      # remove trailing \r\n
        my $pos = "$l:$."; # position
        my @params = split /\t/;
        if ($#params < 4) {
            fail("$pos - invalid line");
            next;
        }
        my ($table, $fq1, $fq1md5, $fq2, $fq2md5, $adpt1, $adpt1md5, $adpt2, $adpt2md5) = @params;
        $table = "fq_$table";
        $table =~ s/\W/_/g; # replace non support char to underline
        $fq1 = "$p/$fq1";
        $fq2 = "$p/$fq2";
        $adpt1 = $adpt1 ? "$p/$adpt1" : "";
        $adpt2 = $adpt2 ? "$p/$adpt2" : "";
        next unless check_md5($fq1, $fq1md5);
        next unless check_md5($fq2, $fq2md5);
        next unless ($adpt1 and $adpt1md5) and check_md5($adpt1, $adpt1md5);
        next unless ($adpt2 and $adpt2md5) and check_md5($adpt2, $adpt2md5);
        # keep _1 the comes first
        ($fq1, $fq2) = ($fq2, $fq1) if ($fq1 gt $fq2);
        upload($table, $fq1, $fq2, $pos, $adpt1, $adpt2);
    }
    close $f;
    debug("<< $l: done");
}

$q->enqueue(undef) for 1..$thread_num;
$_->join() for @runners;

my $l = scalar(@failures);
if ($l > 0) {
    debug("== $l failures:");
    print "$_\n" for (@failures);
    exit 1;
} else {
    debug("== All done.");
    exit 0;
}
