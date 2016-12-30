#!/usr/bin/env perl

use strict;
use warnings;

my $project = shift @ARGV;
my $instance = shift @ARGV;
my $taskname = shift @ARGV;
my $addr = shift @ARGV;
my $token = shift @ARGV;
my $thread = shift @ARGV;

my $cmd = "./download_summary.sh $project $instance $taskname | ./get_logid_list.pl 1 | ./extract_log_info.pl \"<< Tmp file size \" $addr $project $instance $token $thread | sort -k1n |";

open OUT, $cmd;
while (<OUT>) {
    print $_;
}
