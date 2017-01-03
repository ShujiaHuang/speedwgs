#!/bin/bash

if [ "$#" -ne "5" ];then
    echo "Usage: `basename $0` input_table output_table reducer_num split_size <disk|mem>"
    exit 1
fi

input_table=$1
output_table=$2
reducer_num=$3
split_size=$4
use_disk=$5

work_dir=$(cd `dirname $0`; pwd)

job_cmd="jar -resources gene_bin.tar.gz,human_g1k_v37.tar.gz,dbsnp_138.tar.gz,1000G_phase1.indels.tar.gz,Mills_and_1000G_gold_standard.indels.tar.gz com.aliyun.odps.mapred.bridge.streaming.StreamJob -numReduceTasks $reducer_num -partitioner com.aliyun.odps.mapred.lib.KeyFieldBasedPartitioner -jobconf odps.stage.mapper.split.size=$split_size -jobconf stream.map.input.field.separator=\n -jobconf stream.map.output.field.separator=\001 -jobconf map.output.key.field.separator=\001 -jobconf stream.num.map.output.key.fields=3 -jobconf num.key.fields.for.partition=1 -jobconf stream.map.output.key.options=1n,2n,3n -jobconf stream.reduce.input.field.separator=\001 -jobconf odps.stage.mapper.mem=7500 -jobconf odps.stage.reducer.mem=3000 -input $input_table -output $output_table -mapper \"mapper_custompart.pl part.tsv $use_disk\" -reducer \"reducer_gvcf.pl $use_disk\" -file $work_dir/bin/mapper_custompart.pl -file $work_dir/bin/reducer_gvcf.pl -file $work_dir/bin/part.tsv"

odpscmd --config=conf/odps_config.ini -e "create table if not exists $output_table (s string)"
odpscmd --config=conf/odps_config.ini -e "$job_cmd"
