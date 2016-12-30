#!/usr/bin/env perl

use strict;
use warnings;
use Gene_util;

# args
my $part_filename = shift @ARGV;
my $useDisk = shift @ARGV;
die "usage: $0 partfile <disk|mem> [group_header]\n"
	unless defined($useDisk) and defined($part_filename);
my $group_header = shift @ARGV;
unless (defined($group_header)) {
    $group_header = '@RG\tID:HiseqEAAAGAAA-98\tPL:illumina\tPU:150430_I00137_FCH2JTWBBXX_L4_HiseqEAAAGAAA-98\tLB:HiseqEAAAGAAA-98\tSM:YH250\tCN:BGI';
    print STDERR "WARNING: group header not specified, using default:\n$group_header\n";
}

# TODO read $n $q $qs from ARGV
my $n = 0.05;
my $q = 0.5;
my $qs = 33; # 64 or 33

use File::Basename;
my $bin = "gene_bin.tar.gz";

my $tmp_root = -d "/var/run" ? "/var/run" : ".";
if ($useDisk eq "disk") {
    $tmp_root = ".";
}
my $tmp1 = "$tmp_root/$$-fastq-1.fastq";
my $tmp2 = "$tmp_root/$$-fastq-2.fastq";
my $code = 0; # subprocess exit code

# read part file
print STDERR "reading partfile $part_filename ...\n";
my @partinfo; # rname begin end
open PART, "<", $part_filename or die "failed to open $part_filename.\n";
while (<PART>) {
    next if  m/^\s*$/; # ignore empty line
    next if m/^#/;    # ignore comment
    my @item = split /\s+/;
    print STDERR "WARNING: malformed partfile line: $_" unless ($#item == 2);
    push @partinfo, \@item;
}
close PART;
# sort part file
my %idx; # rname -> seq
my $i = 0;
for (1 .. 22, "X", "Y", "MT") {
    $idx{$_} = ++$i;
}
@partinfo = sort { $idx{$a->[0]} <=> $idx{$b->[0]} or $a->[1] <=> $b->[1] } @partinfo;
my $reducer_num = $#partinfo + 1;
print STDERR "$reducer_num reducers according to partfile\n";

sub qc_filter {
    my $qual = shift;
    chomp $qual;
    my $len = length $qual;
    my $qv = 0;
    my $low_q_num = 0;
    my $n_num = 0;
    my $c;
    for my $i (0 .. $len - 1) {
        $c = substr($qual, $i, 1); # current quality char
        $qv = ord($c) - $qs;       # quality value
        $low_q_num++ if ($qv <= 5);
        $n_num++ if ($c eq 'N');
    }
    return $low_q_num / $len > $q || $n_num / $len > $n;
}

# write input stream to 2 local files
open OUT1, ">", $tmp1;
open OUT2, ">", $tmp2;
my $c = 0;
my @f;
while (!eof(STDIN)) {
    for my $i (0 .. 7) {
        $f[$i] = <>;
    }
    if (qc_filter($f[3]) or qc_filter($f[7])) {
        next;
    }
    print OUT1 $f[0], $f[1], $f[2], $f[3];
    print OUT2 $f[4], $f[5], $f[6], $f[7];
    print STDERR "$c records consumed.\n" if (++$c % 10000 == 0);
}
close OUT1;
close OUT2;
print STDERR "done preparing input, $c records in total.\n";

use MIME::Base64;
use List::Util;
use POSIX;

sub get_part_id {
    my $rname = shift;
    my $pos = shift;
    my $id = shift;
    my $step = 0;
    # first shot, and decide direction if miss
    if ($rname eq $partinfo[$id][0]) {
        if ($pos >= $partinfo[$id][1] and $pos <= $partinfo[$id][2]) {
            return $id;
        } elsif ($pos < $partinfo[$id][1]) {
            $step = -1;
        } elsif ($pos > $partinfo[$id][2]) {
            $step = 1;
        } else {
            return -1;
        }
    } elsif ($idx{$rname} < $idx{$partinfo[$id][0]}) {
        $step = -1;
    } elsif ($idx{$rname} > $idx{$partinfo[$id][0]}) {
        $step = 1;
    }
    # pick the direction and go further
    while (1) {
        $id += $step;
        return -1 if ($id < 0 or $id > $#partinfo);
        if ($rname eq $partinfo[$id][0]) {
            if ($pos >= $partinfo[$id][1] and $pos <= $partinfo[$id][2]) {
                return $id;
            } elsif (($step == -1 and $pos > $partinfo[$id][2]) or
                     ($step == 1 and $pos < $partinfo[$id][1])) {
                return -1;
            }
        }
    }
}

my @std_ref= Gene_util::get_std_ref();

my %chr_seq; # rname -> sequence(start from 1)
my $sum = 0;
my %index; # rname -> offset
for (my $i = 0; $i <= $#std_ref ; $i += 2) {
    $index{$std_ref[$i]} = $sum;
    $sum += $std_ref[$i + 1];
    $chr_seq{$std_ref[$i]} = int($i / 2) + 1;
}
my $each_part = ceil($sum / $reducer_num);

my $delim = "\x01";
my $header = "";
my %reducer_header;

my $pid = fork();

if ($pid == 0) { # child process do the computing
    open BWA, "$bin/bwa mem -M -a -R '$group_header' human_g1k_v37.tar.gz/human_g1k_v37.fasta $tmp1 $tmp2 |";
    # ERR194147.1     97      MT      16503   60      67M34S  =       1       -16403  GGTTCCTACTTCAGGGTCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACG   CC@FFFFFHHHHHJJJFHIIJJJJJJIHJIIJJJJJJJJIIGIJJIJJJIJJJIJIJJJJJJJJJJIJHHHHFFFDEEEEEEEEDDDCDDEEDDDDDDDDD   NM:i:0  MD:Z:67 AS:i:67 XS:i:33 RG:Z:HiseqEAAAGAAA-98   SA:Z:17,22020694,+,34S67M,34,4;
    while (<BWA>) {
        if (m/^\@/) {
            $header .= $_;
            next;
        }
        my (undef, undef, $rname, $pos, undef) = split(/\t/, $_, 5);
        next unless (exists $index{$rname});                    # skip rname = GLxxx|*
        my $part = floor(($index{$rname} + $pos) / $each_part); # use even split as start part id
        $part = $part > $reducer_num - 1 ? $reducer_num - 1 : $part;
        $part = get_part_id($rname, $pos, $part);               # search for proper part id
        next if ($part == -1);
        unless (exists $reducer_header{$part}) {
            $reducer_header{$part} = 1;
            # use chr_seq = 0 (smallest seq), pos = 0 as sign of header
            print join($delim, $part, 0, 0, encode_base64($header, "")), "\n";
        }
        print join($delim, $part, $chr_seq{$rname}, $pos, $_);
    }
    close BWA;
    $code = $? >> 8;

    unlink $tmp1, $tmp2; # clean tmp files
    exit($code);
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
