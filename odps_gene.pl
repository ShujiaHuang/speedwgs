#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;

my $output_field = "output_content"; # field name of tmp tables
my $volume_name = "bgi";
my $work_dir = cwd();

my $path= $ENV{'PATH'};
$ENV{'PATH'} = "$path:$work_dir/odpscmd/bin:$work_dir/tattoo/bin";

sub print_usage {
print STDERR "Usage:
    -setup resource_dir
    -upload-fastq file1 file2 table_name conf_file thread_num
    -run-sql \"sql-command\"
    -haplotype target_table input_table reducer_number group_header
    -union target_table input_table1,input_table2,...
    -genotype target_table input_table reducer_number
    -combine-vcf target_table input_table
";
}

if ($#ARGV == -1) {
    &print_usage();
    exit;
}

sub read_arg {
    my $param = shift;
    my $arg = shift @ARGV; # read arg from ARGV
    die "$param not specified!\n" unless defined $arg;
    return $arg;
}

my $cmd = lc(&read_arg("command"));

if ($cmd eq "-run-sql") {
    &run_sql();
} elsif ($cmd eq "-haplotype") {
    &run_getvcf_haplotype();
} elsif ($cmd eq "-union") {
    &run_union_gvcfs();
} elsif ($cmd eq "-genotype") {
    &run_genotype_gvcfs();
} elsif ($cmd eq "-combine-vcf") {
    &run_combine_vcf();
} elsif ($cmd eq "-upload-fastq") {
    &run_upload_fastq();
} elsif ($cmd eq "-setup") {
    &run_setup();
} else {
    &print_usage();
    exit;
}

sub run_sql {
    my $sql = &read_arg("sql");
    &run_odps_cmd($sql);
}

sub run_getvcf_haplotype {
    my $target_table = &read_arg("target table");
    my $input_table = &read_arg("input table");
    my $reducer_num = &read_arg("reducer number");
    my $group_header = &read_arg("group header");
    &run_odps_cmd("CREATE TABLE $target_table ($output_field STRING)");
    my $haplotype_cmd = &get_job_cmd(
        $input_table,
        $target_table,
        256,
        $reducer_num,
        11264,
        6144,
        "mapper_custompart.pl part.tsv disk $group_header",
        "reducer_gvcf.pl disk gatk3.5 gvcf $input_table",
        3,
        "1n,2n,3n",
        "\\n",
        "part.tsv"
        );
    &run_odps_cmd($haplotype_cmd);
}

sub get_union_cmd {
    my ($target_table, @input_tables) = @_;
    my $cmd = "CREATE TABLE $target_table AS SELECT table_name, $output_field FROM (";
    my $union_str = "";
    for my $gvcf_table (@input_tables) {
        $cmd .= "$union_str SELECT \\\"$gvcf_table\\\" AS table_name, $output_field from $gvcf_table";
        $union_str = " UNION ALL";
    }
    $cmd .= ") t";
    return $cmd;
}

sub run_union_gvcfs {
    my $target_table = &read_arg("target table");
    my @input_tables = split /,/, &read_arg("input tables");
    my $union_cmd = &get_union_cmd($target_table, @input_tables);
    &run_odps_cmd($union_cmd);
}

sub run_genotype_gvcfs {
    my $target_table = &read_arg("target table");
    my $input_table = &read_arg("input table");
    my $reducer_num = &read_arg("reducer number");
    &run_odps_cmd("CREATE TABLE $target_table ($output_field STRING)");
    my $genotype_cmd = &get_job_cmd(
        $input_table,
        $target_table,
        256,
        $reducer_num,
        3072,
        6144,
        "mapper_genotype.pl",
        "reducer_genotype.pl",
        4,
        "1n,2,3n,4n",
        "\\001"
        );
    &run_odps_cmd($genotype_cmd);
}

sub run_combine_vcf {
    my $target_table = &read_arg("target table");
    my $input_table = &read_arg("input table");
    &run_odps_cmd("CREATE TABLE $target_table ($output_field STRING)");
    my $combine_cmd = &get_job_cmd(
        $input_table,
        $target_table,
        256,
        1,
        3072,
        6144,
        "mapper_combinevcf.pl",
        "reducer_combinevcf.pl",
        3,
        "1n,2n,3n",
        "\\001"
        );
    &run_odps_cmd($combine_cmd);
}

sub run_upload_fastq {
    my $file1 = &read_arg("file1");
    my $file2 = &read_arg("file2");
    my $table_name = &read_arg("table name");
    my $conf_file = &read_arg("odps conf file");
    my $thread_num = &read_arg("thread number");
    my $cmd = "odpscmd --config=$conf_file -e 'create table if not exists $table_name (id1 string, seq1 string, opt_id1 string, quality1 string, id2 string, seq2 string, opt_id2 string, quality2 string)'";
    0 == system($cmd) or die "failed to create table $table_name\n";
    $cmd = "$work_dir/data_upload/merge_fastq.pl $file1 $file2 $table_name $conf_file $thread_num";
    0 == system($cmd) or die "command failed: $cmd\n";
}

sub run_setup {
    my $resource_dir = &read_arg("resource file dir");

    my $bin_cmd = "add archive $work_dir/bin/bin.tar.gz as gene_bin.tar.gz";
    &run_odps_cmd($bin_cmd);
    my $mkv_cmd = "fs -mkv $volume_name";
    &run_odps_cmd($mkv_cmd);

    my @ref_resources = ("human_g1k_v37.tar.gz",
                         "dbsnp_138.tar.gz",
                         "1000G_phase1.indels.tar.gz",
                         "Mills_and_1000G_gold_standard.indels.tar.gz");
    for my $ref_archive (@ref_resources) {
        &add_volume_archive($resource_dir, $ref_archive);
    }
}

sub add_volume_archive {
    my $dir = shift;
    my $file = shift;
    my $abs_path = "$dir/$file";
    my $volume_part = substr($file, 0, index($file, "."));
    print STDERR "Uploading resource file: $abs_path, this may take a while...\n";
    my $cmd = "fs -put $abs_path /$volume_name/$volume_part; add volumearchive /$volume_name/$volume_part/$file as $file";
    &run_odps_cmd($cmd);
    print STDERR "Uploade success! File: $abs_path\n";
}

sub run_odps_cmd {  # usage: $ret = run_odps_cmd(command)
    my $cmd = shift;
    0 == system("odpscmd --config=conf/odps_config.ini -e \"$cmd\"") or die "command failed: $cmd\n";
}

sub get_job_cmd {
    my $input_table = shift;
    my $target_table = shift;
    my $map_split_size = shift;
    my $reducer_num = shift;
    my $map_memory = shift;
    my $reduce_memory = shift;
    my $mapper = shift;
    my $reducer = shift;
    my $map_output_key_fields = shift;
    my $map_output_key_options = shift;
    my $map_input_seperator = shift;
    my $aux_file = shift; # optional

    my $append_file = defined($aux_file) ? "-file $work_dir/bin/$aux_file" : "";
    my ($mapper_script) = grep length, split /\s+/, $mapper;
    my ($reducer_script) = grep length, split /\s+/, $reducer;
    return "
jar -resources gene_bin.tar.gz,human_g1k_v37.tar.gz,dbsnp_138.tar.gz,1000G_phase1.indels.tar.gz,Mills_and_1000G_gold_standard.indels.tar.gz
    com.aliyun.odps.mapred.bridge.streaming.StreamJob
    -numReduceTasks $reducer_num
    -partitioner com.aliyun.odps.mapred.lib.KeyFieldBasedPartitioner
    -jobconf odps.stage.mapper.split.size=$map_split_size
    -jobconf stream.map.input.field.separator=$map_input_seperator
    -jobconf stream.map.output.field.separator=\\001
    -jobconf map.output.key.field.separator=\\001
    -jobconf stream.num.map.output.key.fields=$map_output_key_fields
    -jobconf num.key.fields.for.partition=1
    -jobconf stream.map.output.key.options=$map_output_key_options
    -jobconf stream.reduce.input.field.separator=\\001
    -jobconf odps.stage.mapper.mem=$map_memory
    -jobconf odps.stage.reducer.mem=$reduce_memory
    -input $input_table
    -output $target_table
    -mapper \\\"$mapper\\\"
    -reducer \\\"$reducer\\\"
    -file $work_dir/bin/$mapper_script
    -file $work_dir/bin/$reducer_script
    -file $work_dir/bin/Gene_util.pm
    $append_file";
}
