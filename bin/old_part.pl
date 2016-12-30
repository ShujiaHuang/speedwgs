#!/usr/bin/env perl

use strict;
use warnings;

use POSIX;

my $reducer_num = shift @ARGV;
die "usage: $0 reducer_num\n" unless defined($reducer_num);

my @seq = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
           "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT");

my %std_ref = (
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
    "MT" => 16569);

my $sum = 0;
for my $c (@seq) {
    $sum += $std_ref{$c};
}
my $each_part = ceil($sum / $reducer_num);
print ">> $each_part for each reducer, total length = $sum\n";

my $id = 0;
my $begin = 0;
my $last = "";
for my $c (@seq) {
    my $len = $std_ref{$c};
    while ($begin + $each_part < $len) {
        print "$id\t", $begin < 0 ? "$last,$c:1-" : "$c:$begin-", $begin + $each_part - 1, "\n";
        $begin += $each_part;
        $id++;
    }
    $last = "$c:$begin-$std_ref{$c}";
    $begin = $begin - $std_ref{$c};
}
# TODO last reducer
# my $c = $seq[-1];
# print "$id\t", $begin < 0 ? "$last,$c:1-" : "$c:$begin-", $begin + $each_part, "\n";
