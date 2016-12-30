#!/usr/bin/env perl
use strict;
use warnings;
use POSIX qw(mkfifo);

my $f1 = shift @ARGV;
my $f2 = shift @ARGV;
my $table_name = shift @ARGV;
my $conf_file = shift @ARGV;
my $thread_num = shift @ARGV;
die "usage: $0 file1 file2 table_name conf_file thread_num [adapter...]\n"
    unless(defined($f1) and defined($f2) and defined($table_name) and defined($conf_file) and defined($thread_num));
my %adpt;
while (my $f = shift @ARGV) {
    next unless $f;
    open ADPT, "<", $f or die "failed to open adapter file $f\n";
    while (<ADPT>) {
        next if /^#/;
        my $id = (split /\s+/)[0];
        $id =~ s/\/[12]$//;
        $adpt{$id}++;
    }
    close ADPT;
}
sub in_adpt {
    my $id = shift;
    $id = (split /\s+/, $id)[0];
    $id =~ s/^@//;
    $id =~ s/\/[12]$//;
    return exists $adpt{$id};
}

my $tmp_file = "$$-tmp";
mkfifo($tmp_file, 0700);

sub get_in {
    my $file = shift;
    if ($file =~ m/gz$/) {
        return "zcat $file |";
    }
    return "cat $file |";
}

my $pid = fork();
if ($pid == 0) {
    open IN1, &get_in($f1);
    open IN2, &get_in($f2);
    open OUT, ">", $tmp_file;
    my (@l, @r);
    while (!eof(IN1) and !eof(IN2)) {
        $l[$_] = <IN1> for (0 .. 3);
        $r[$_] = <IN2> for (0 .. 3);
        next if (in_adpt($l[0]) or in_adpt($r[0]));
        print OUT @l, @r;
    }
    close IN1;
    close IN2;
    close OUT;
    exit;
} elsif ($pid) {
    my $ret = system("tattoo -w $thread_num --config $conf_file -t $table_name --quiet read-text $tmp_file --num-lines 8");
    if ($ret == 0) {
        wait();
    } else {
        kill 'KILL', $pid;
    }
}

unlink($tmp_file);
