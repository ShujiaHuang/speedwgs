#!/usr/bin/env perl

use strict;
use warnings;
use MIME::Base64;
use POSIX qw(WNOHANG);
use POSIX ":sys_wait_h";
use Gene_util;

Gene_util::check_eof();

my $bin = "gene_bin.tar.gz";
my $delim = "\x01";

my @vcf_files;

my $fh;
my $idx_fh;
while (<>) { # prepare input vcf files
    my ($vcf_id, $chr_seq, $pos, $vcf_content) = split /$delim/;
    my $file_name = "$vcf_id.vcf";
    if (!$vcf_files[$#vcf_files] or
        ($vcf_files[$#vcf_files] ne $file_name)) {
        if ($vcf_files[$#vcf_files]) {
            close $fh;
        }
        open ($fh, '>', $file_name);
        push @vcf_files, $file_name;
    }
    if ($chr_seq == -1 and $pos == -1) {
        print $fh decode_base64($vcf_content);
        next;
    }
    if ($chr_seq == -1 and $pos == 0) {
        my $idx_file_name = "$file_name.idx";
        open $idx_fh, '>', $idx_file_name;
        print $idx_fh decode_base64($vcf_content);
        close $idx_fh;
        next;
    }
    print $fh $vcf_content;
}
close $fh;

sub get_combine_cmd {
    my $cmd = "$bin/gatk3.5 -T CombineVariants -R human_g1k_v37.tar.gz/human_g1k_v37.fasta ";
    for my $file (@_) {
        $cmd .= "--variant $file ";
    }
    return $cmd." --genotypemergeoption UNSORTED |";
}

my $pid = fork();
if ($pid == 0) {
    my $code = 0;
    my %chr_seq = Gene_util::get_chr_seq();
    my $cmd = get_combine_cmd(@vcf_files);
    my $vcf_header = "";
    my $out_header_printed = 0;
    open OUT, $cmd;
    while (<OUT>) {
        if (m/^#/) {
            $vcf_header .= $_;
            next;
        }
        if (!$out_header_printed) {
            print join($delim, -1, -1, encode_base64($vcf_header, "")), "\n";
            $out_header_printed = 1;
        }
        my ($rname, $pos, $others) = split("\t", $_, 3);
        print join($delim, $chr_seq{$rname}, $pos, $_);
    }
    close OUT;
    $code = $? >> 8;
    print STDERR "<< Tmp file size ";
    system("du -s 1>&2");
    for my $file (@vcf_files) {
        unlink($file, "$file.idx");
    }
    exit($code);
} elsif ($pid) {
    my $c = 0;
    while (1) {
        my $child = waitpid($pid, WNOHANG);
        if ($child != 0) {
            $? == 0 ? print STDERR "everything seems good, exit\n" : print STDERR "something bad happened.\n";
            exit($? >> 8);
        }
        if ($c++ > 60 * 5) { # default heartbreak timeout is 10 min
            print STDERR "heartbeat\n";
            $c = 0;
        }
        sleep(1);
    }
}
