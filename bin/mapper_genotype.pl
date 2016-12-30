#!/usr/bin/env perl

use strict;
use warnings;

my $delim = "\x01";
while (<>) {
    my ($table_name, $part_id, $others) = split /$delim/, $_, 3;
    print join($delim, $part_id, $table_name, $others); # let part_id be the first field
}
