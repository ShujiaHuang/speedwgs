#!/usr/bin/env perl

use strict;
use warnings;
use threads;
use Thread::Queue;
use HTTP::Tiny;

my $grep_hint = shift @ARGV;
my $ip_addr = shift @ARGV;
my $project = shift @ARGV;
my $instance_id = shift @ARGV;
my $token = shift @ARGV;
my $thread_num = shift @ARGV;

die "Usage: $0 grep_hint ip_address project instance_id authoration_token\n"
    unless defined($grep_hint) and defined($ip_addr) and defined($project) and defined($instance_id) and defined($token);

$thread_num = 1 unless defined($thread_num); # default 1 thread

my $begin_hint = ">> start worker: ";
my $end_hint = "everything seems good, exit";

my @runner_threads;
my $q = Thread::Queue->new();
my $ht = HTTP::Tiny->new();

for (1..$thread_num) {
    push @runner_threads, async {
        while (my $input = $q->dequeue()) {
            my ($worker_idx, $log_id) = split /,/, $input;
            chomp($log_id);
            my $log_url = "http://$ip_addr/api/projects/$project/instances/$instance_id?log&logtype=Stderr&size=5000000&id=$log_id&authorization_token=$token";
            my $response = $ht->get($log_url);
            unless ($response->{success}) { # if get nothing, print err & next
                print STDERR "Get log failed. Worker index: $worker_idx, Url: $log_url\n";
                next;
            }
            my @lines = split /\n/, $response->{content};
            print $worker_idx, " ", join(" ", &process_log($worker_idx, $grep_hint, @lines)), "\n";
        }
    }
}

while (<>) {
    $q->enqueue($_);
}
$q->enqueue(undef) for 1..$thread_num;
$_->join() for @runner_threads;

sub process_log {
    my $worker_idx = shift;
    my $grep_hint = shift;
    my @lines = @_;
    my $is_begin = 0;

    my @collector;
    for my $line (@lines) {
        unless ($is_begin) {
            my $start = index($line, $begin_hint);
            if ($start != -1 and substr($line, $start + length($begin_hint)) == $worker_idx) { # check begin hint
                $is_begin = 1;
            }
            next;
        }
        last if($line =~ /$end_hint/); # break if end hint

        my $start_idx = index($line, $grep_hint);
        next unless $start_idx != -1;
        $start_idx += length($grep_hint);
        my $end_idx = index($line, " ", $start_idx);
        push @collector, substr($line, $start_idx, $end_idx);
    }
    return @collector;
}
