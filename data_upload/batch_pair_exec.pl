#!/usr/bin/env perl
use strict;
use warnings;
use threads;
use Thread::Queue;

my %files;
my $thread_num = shift @ARGV;
my $exec_script = shift @ARGV;

die "usage: $0 thread_num exec_script [exec_script args]\n"
    unless (defined($exec_script) and defined($thread_num));

sub get_pair {
    my $key = shift;
    my $feature = "_R";
    my $i = rindex($key, $feature) + 2;
    if (substr($key, $i, 1) eq "1") {
        substr($key, $i, 1) = "2"; # switch key to pair->second
        return $key; # only return pair_second if this is pair_first
    }
    return; # return nothing if this is not pair_first
}

while (my $line = <STDIN>) {
    chomp($line);
    $files{$line} = 1;
}

my @runner_threads;
my $q = Thread::Queue->new();

for (1..$thread_num) {
    push @runner_threads, async {
        while (my $pair_first = $q->dequeue()) {
            my $pair_second = &get_pair($pair_first);
            my $cmd = "./$exec_script $pair_first $pair_second @ARGV";
            0 == system($cmd) or die "$cmd\n";
        }
    };
}

for my $k (keys %files) {
    my $pair = &get_pair($k);
    if ($pair and $files{$pair}) { # if pair_second exists, exec pair script
        $q->enqueue($k);
    }
}

$q->enqueue(undef) for 1..$thread_num;
$_->join() for @runner_threads;
