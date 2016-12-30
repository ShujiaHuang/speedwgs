#!/usr/bin/env perl

use strict;
use warnings;
use MIME::Base64;
use POSIX qw(WNOHANG);
use POSIX;
use POSIX ":sys_wait_h";
use Gene_util;

Gene_util::check_eof(); # Exit if no input

my $bin = "gene_bin.tar.gz";
my $delim = "\x01";

my %header_printed; # table_name -> bool
my @gvcf_files;
my $vcf_id;
my $fh;
while (<STDIN>) {
    my ($part_id, $table_name, $chr_seq, $pos, $gvcf) = split /$delim/;
    $vcf_id = $part_id;
    # The last of gvcf_files is being written, filename: $table_name.g.vcf
    my $file_name = "$table_name.g.vcf";
    if (!$gvcf_files[$#gvcf_files] or
        ($gvcf_files[$#gvcf_files] ne $file_name)) { # A new table comes

        if ($gvcf_files[$#gvcf_files]) { # close the old one if it exists
            close $fh;
        }
        open ($fh, '>', $file_name);
        push @gvcf_files, $file_name;
    }
    if ($chr_seq == -1 and $pos == -1) {
        if (!$header_printed{$table_name}) {
            print $fh decode_base64($gvcf);
            $header_printed{$table_name} = 1;
        }
    } else {
        print $fh $gvcf;
    }
}
close $fh;

print STDERR ">> start worker: $vcf_id\n";

sub get_genotype_cmd {
    my $tmp_vcf = shift;
    my $cmd = "$bin/gatk3.5 -T GenotypeGVCFs -R human_g1k_v37.tar.gz/human_g1k_v37.fasta ";
    for my $gvcf (@_) {
        $cmd .= "--variant $gvcf ";
    }
    return $cmd."-o $tmp_vcf";
}

my $pid = fork();

if ($pid == 0) { # child process do the computing
    my %chr_seq = Gene_util::get_chr_seq();

    my $tmp_vcf_out = "tmp.vcf";
    my $cmd = get_genotype_cmd($tmp_vcf_out, @gvcf_files);
    Gene_util::system2($cmd);

    open VCF_OUT, '<', $tmp_vcf_out;
    my $out_header = "";
    my $out_header_printed = 0;
    while (<VCF_OUT>) {
        if (m/^#/) {
            $out_header .= $_;
            next;
        }
        if (!$out_header_printed) {
            print join($delim, $vcf_id, -1, -1, encode_base64($out_header, "")), "\n";
            $out_header_printed = 1;
        }
        my ($rname, $pos, $others) = split("\t", $_, 3);
        print join($delim, $vcf_id, $chr_seq{$rname}, $pos, $_);
    }
    close VCF_OUT;

    if ($out_header_printed) { # if header is printed, also print idx file
        open IDX_OUT, '<', "$tmp_vcf_out.idx"; # gatk will generate this file automatically
        my $idx_content = "";
        while (<IDX_OUT>) {
            $idx_content .= $_;
        }
        close IDX_OUT;
        print join($delim, $vcf_id, -1, 0, encode_base64($idx_content, "")), "\n";
    }

    print STDERR "<< Tmp file size ";
    system("du -s 1>&2");
    for my $file (@gvcf_files) {
        unlink($file, "$file.idx");
    }
    unlink($tmp_vcf_out, "$tmp_vcf_out.idx");
    exit(0);
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
