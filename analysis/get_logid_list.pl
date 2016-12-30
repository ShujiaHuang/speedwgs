#!/usr/bin/env perl

use strict;
use warnings;
use JSON;
use Data::Dumper;

my $stage_idx = shift @ARGV;

die "Usage: $0 stage_idx\n"
    unless defined($stage_idx);

my $json_str = "";
my $start = 0;
while (<>) {
    if (m/^{/) {
        $start = 1;
    }
    if ($start) {
        $json_str .= $_;
    }
}

sub get_worker_idx {
    my $worker_id = shift;
    my $start = rindex($worker_id, "#") + 1;
    my $end = rindex($worker_id, "_");
    return substr($worker_id, $start, $end - $start);
}

my $instance_log = decode_json($json_str);
my @workers = @{$instance_log->{"Instance"}{"Stages"}[$stage_idx]{"Workers"}};
for my $worker (@workers) {
    print &get_worker_idx($worker->{"WorkerID"}) . "," . $worker->{"LogID"},"\n";
}
