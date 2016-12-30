#!/usr/bin/env perl

use strict;
use warnings;

while(<>) {
    if (m/^ALIYUN/) {
        my ($a, $b) = split /:/;
        system("./odps_gene.pl -run-sql 'drop table $b'");
    }
}
