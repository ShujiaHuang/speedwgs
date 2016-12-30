#!/usr/bin/env perl

use strict;
use warnings;

use MIME::Base64;

my $table = shift @ARGV;
my $output = shift @ARGV;
unless (defined($table) and defined($output)) {
    print "usage: $0 table file [odps_config.ini]\n";
    print "\tthis script will download vcf table from ODPS to local file.\n";
    exit 1;
}
my $odps_config = shift @ARGV;
my $odps_param = defined($odps_config) ? "--config=$odps_config" : "";

my $tmp = "$$-tmp.vcf";
use POSIX qw(mkfifo);
mkfifo($tmp, 0700) or die "failed to mkfifo $tmp\n";
my $tunnel_pid = fork();
die "FATAL: failed to fork process for downloading table\n" unless defined($tunnel_pid);
if ($tunnel_pid == 0) {
    print STDERR ">> downloading table $table info $tmp ...\n";
    0 == system("odpscmd/bin/odpscmd $odps_param -e 'tunnel d $table $tmp'") or die "failed to download $table.\n";
    #0 == system("odpscmd $odps_param -e 'tunnel d $table $tmp'") or die "failed to download $table.\n";
    print STDERR "<< download complete\n";
    exit 0;
}
my $sort_pid = fork();
die "FATAL: failed to fork process for sorting result\n" unless defined($sort_pid);
if ($sort_pid == 0)  {
    print STDERR ">> sort and handle header ...\n";
    open VCF, "sort -k1,1n -k2,2n -t \$'\\x01' -T . $tmp |";
    open OUT, ">", $output;
    my $header_printed = 0;
    while (<VCF>) {
        my ($reducer_id, $ln, $content) = split (/\001/);
        if ($reducer_id == -1 and $ln == -1) {
            if (!$header_printed) {
                print OUT decode_base64($content);
                $header_printed = 1;
            }
        } else {
            print OUT $content;
        }
    }
    close OUT;
    close VCF;
    unlink $tmp;
    print STDERR "<< sort complete\n";
}

my $code = 0;
waitpid($tunnel_pid, 0);
$code = $? >> 8;
exit $code unless ($code == 0);
waitpid($sort_pid, 0);
$code = $? >> 8;
exit $code unless ($code == 0);
