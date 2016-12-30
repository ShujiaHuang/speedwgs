#!/usr/bin/env perl
use strict;
use warnings;

my $f1 = shift @ARGV;
my $f2 = shift @ARGV;

die "Usage: $0 pair_first pair_second\n"
    unless (defined($f1) and defined($f2));

sub get_in {
    my $file = shift;
    if ($file =~ m/gz$/) {
        return "zcat $file |";
    }
    return "cat $file |";
}

sub die_error {
    my $content1 = shift;
    my $content2 = shift;
    print STDOUT $f1,"\n";
    print STDOUT $f2,"\n";
    print STDERR "-----\nNot match!\nf1: $f1\nf2: $f2\n$content1\n$content2\n-----\n";
}

open IN1, &get_in($f1);
open IN2, &get_in($f2);

while (!eof(IN1) and !eof(IN2)) {
    my $line_11 = <IN1>;
    my $line_12 = <IN1>;
    my $line_13 = <IN1>;
    my $line_14 = <IN1>;

    my $line_21 = <IN2>;
    my $line_22 = <IN2>;
    my $line_23 = <IN2>;
    my $line_24 = <IN2>;

    if (substr($line_11, 0, 1) ne substr($line_21, 0, 1)) {
        die_error($line_11.$line_12.$line_13.$line_14, $line_21.$line_22.$line_23.$line_24);
        last;
    }
    if ($line_13 ne $line_23) {
        die_error($line_11.$line_12.$line_13.$line_14, $line_21.$line_22.$line_23.$line_24);
        last;
    }
}

close IN1;
close IN2;
