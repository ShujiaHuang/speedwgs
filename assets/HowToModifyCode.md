# 代码阅读指南

odps-gene-compute 这个项目的代码由大量的 ODPS Streaming 作业构成。本文着重讲述如何学习、理解这些代码。

## 准备知识

### ODPS Streaming 作业

ODPS Streaming 作业和 Hadoop Streaming 类似，是一种 MapReduce 作业，允许 Mapper 和 Reducer 使用任意语言完成。

输入表的数据会被转化成指定分隔符分割的字节流，通过 STDIN 传递给 Mapper/Reducer。Mapper/Reducer 处理后的数据按指定分隔符输出到 STDOUT，并被存放到输出表中。

阅读 `odps_productdoc_mr_streaming.pdf` 可以获得更为详细的信息。

### odps_gene.pl

启动 ODPS Streaming 作业可以在 odpscmd 中完成，但是命令行很长，类似这样：

```
jar -resources gene_bin.tar.gz,human_g1k_v37.tar.gz,dbsnp_138.tar.gz,1000G_phase1.indels.tar.gz,Mills_and_1000G_gold_standard.indels.tar.gz
    com.aliyun.odps.mapred.bridge.streaming.StreamJob
    -numReduceTasks 899
    -partitioner com.aliyun.odps.mapred.lib.KeyFieldBasedPartitioner
    -jobconf odps.stage.mapper.split.size=256
    -jobconf stream.map.input.field.separator=\n
    -jobconf stream.map.output.field.separator=\001
    -jobconf map.output.key.field.separator=\001
    -jobconf stream.num.map.output.key.fields=3
    -jobconf num.key.fields.for.partition=1
    -jobconf stream.map.output.key.options=1n,2n,3n
    -jobconf stream.reduce.input.field.separator=\001
    -jobconf odps.stage.mapper.mem=8000
    -jobconf odps.stage.reducer.mem=6500
    -input hdjy_na12878_fastq_merge
    -output tmp_gvcf_hdjy_na12878_fastq_merge_1473316239
    -mapper "mapper_custompart.pl part.tsv disk '@RG\tID:HiseqEAAAGAAA-98\tPL:illumina\tPU:150430_I00137_FCH2JTWBBXX_L4_HiseqEAAAGAAA-98\tLB:HiseqEAAAGAAA-98\tSM:hdjy_na12878_fastq_merge\tCN:BGI'"
    -reducer "reducer_gvcf.pl disk gatk3.5 gvcf hdjy_na12878_fastq_merge"
    -file /apsarapangu/disk1/tianli.tl/odps-gene-compute2/bin/mapper_custompart.pl
    -file /apsarapangu/disk1/tianli.tl/odps-gene-compute2/bin/reducer_gvcf.pl
    -file /apsarapangu/disk1/tianli.tl/odps-gene-compute2/bin/Gene_util.pm
    -file /apsarapangu/disk1/tianli.tl/odps-gene-compute2/bin/part.tsv
```

这样的命令在 odpscmd 中靠手工正确输入是比较困难的。因此我们封装了 odps_gene.pl 这个脚本，使启动作业尽量简化。

## 作业解析: mapper_custompart.pl + reducer_gvcf.pl

mapper_custompart.pl

- Line 73-91: 首先，从 STDIN 读入表数据，写到两个本地文件，并同时做 qc filter
- Line 152-176: 调用 bwa，对于 bwa 输出的每一行，提取染色体(rname)、位置(pos)两个字段，根据这两个字段计算应该分配给哪个 Reducer(part)，并将 part，rname，pos 以及原输出行写到 STDOUT。这里有个比较 tricky 的地方，就是对 header 的处理。header 是多行字符串，这里先把所有 header 拼接到一个字符串里，然后做 base64 之后输出，其 rname 和 pos 都是 0 以保证在 Reducer 端输入的时候 header 排序位于最前端。

Mapper 的输出会被发送给 Reducer。ODPS 根据启动命令的参数将输出分发给不同 Reducer，并排序。
```
    -jobconf stream.map.output.field.separator=\001    # Mapper 输出的字段以 \x01 分隔
    -jobconf map.output.key.field.separator=\001       # Mapper 输出的字段中用于分发和排序的 key 也使用 \x01 分隔
    -jobconf stream.num.map.output.key.fields=3        # 输出的 key 字段有 3 个（part, rname, pos）
    -jobconf num.key.fields.for.partition=1            # 用于分发的 key 有 1 个（于是剩下的两个用于排序）
    -jobconf stream.map.output.key.options=1n,2n,3n    # 这三个字段都是数值型（否则字符串会被 hash 并分发，数值则直接根据值本身来指定 Reducer）
    -jobconf stream.reduce.input.field.separator=\001  # Reducer 的输入字段以 \x01 分隔
```

reducer_gvcf.pl

- Line 50-92: 从 STDIN 读入数据（包括对 rname=0, pos=0 的头做 base64 decode），并管到给 samtools rmdup，写出去重后的 bam 文件
- Line 118-130: 计算并生成输入数据的 interval
- Line 133-154: 陆续调用 samtools index、gatk IndelRealigner 等
- Line 170-192: 将最终结果（reducer id, rname, pos, 原始 vcf 行）输出到 STDOUT

## 改写作业

ODPS 的存储都是表，如果需要 BAM 落盘，实际上是 SAM 落表。

可以仿照 mapper_custompart.pl 写一个 map-only（没有 reducer）的作业，将 bwa 处理完的数据直接写到输出表。

然后改写现有的 mapper_custompart.pl，将其中的 bwa 逻辑去掉，直接输出出去就可以了。

新的作业写好之后，可以参考 odps_gene.pl 的现有代码（这部分代码还是挺清晰的，就不赘述了），把启动新作业的命令添加进去。

需要注意一点：`mapper_custompart.pl` 和 `reducer_gvcf.pl` 中都进行了 fork 动作，子进程负责具体的运算，父进程负责监控子进程并定时向 STDERR 输出心跳信息。这是因为 ODPS Streaming 作业会对 STDOUT 或 STDERR 监控，如果长时间没有输出（10 分钟），会认为进程陷入死锁而杀掉当前进程并重启新的进程。我们调用的第三方代码无法保证定时有输出，因此代码中采用这种方法来通知 ODPS 进程仍然在健康运行。新写的 Streaming 作业应遵循这一基本框架。
