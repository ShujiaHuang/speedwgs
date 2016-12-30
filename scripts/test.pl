#!/usr/bin/env perl

use strict;
use warnings;
use threads;
use Thread::Queue;

#my $timestamp = time();
my $timestamp = 7;
my $thread_num = 3;
my @runner_threads;

sub get_sample_names {
    my @samples;
    for my $idx (@_) {
        push @samples, "hdjy_na$idx"."_fastq_merge";
    }
    return @samples;
}

sub get_tmp_gvcf_name {
    my $sample_name = shift;
    return "tmp_gvcf_$sample_name"."_$timestamp";
}

sub call_odps_gene {
    my $cmd = shift;
    0 == system($cmd) or die "failed $cmd";
}

############# job 1: for each sample, haplotype #################
my @input_sample_tables = &get_sample_names(12878,12891,12892);
#
#my $q = Thread::Queue->new();
#
#for (1..$thread_num) {
#    push @runner_threads, async {
#        while (my $sample_table = $q->dequeue()) {
#            my $output_table = &get_tmp_gvcf_name($sample_table);
#            &call_odps_gene("./odps_gene.pl -haplotype $output_table $sample_table 899 \\\'\@RG\\\\tID:HiseqEAAAGAAA-98\\\\tPL:illumina\\\\tPU:150430_I00137_FCH2JTWBBXX_L4_HiseqEAAAGAAA-98\\\\tLB:HiseqEAAAGAAA-98\\\\tSM:YH250\\\\tCN:BGI\\\'");
#        }
#    };
#}
#
#
#for my $sample_table (@input_sample_tables) {
#    $q->enqueue($sample_table);
#}
#$q->enqueue(undef) for 1..$thread_num;
#$_->join() for @runner_threads;

########### job 2: union ####################

my @union_input_tables;
for my $sample (@input_sample_tables) {
    push @union_input_tables, &get_tmp_gvcf_name($sample);
}

my $union_output = "tmp_merged_gvcfs_$timestamp";
&call_odps_gene("./odps_gene.pl -union $union_output ".join(',', @union_input_tables));

############ job 3: genotype #############
my $genotype_output = "tmp_genotype_output_$timestamp";
&call_odps_gene("./odps_gene.pl -genotype $genotype_output $union_output 899");

########### job 4: combine ##############
my $combined_vcf_output = "combined_vcf_output_$timestamp";
&call_odps_gene("./odps_gene.pl -combine-vcf $combined_vcf_output $genotype_output");
