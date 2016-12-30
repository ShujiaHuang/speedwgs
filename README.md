# 基于MaxCompute的WGS极速分析流程简要使用说明

极速全基因组数据分析流程名称：**speedwgs**

注意事项：
- 该流程仅实现从原始fastq到输出包含所有候选变异的VCF结果。中间不会输出任何比对（bam/cram）文件。

## 环境准备

### ODPS 工具准备

开通MaxCompute服务之后，通过命令行工具`odpscmd`来对其进行操作。流程包中 odpscmd 目录下已经准备了一份配置好的odpscmd，`odpscmd/bin/odpscmd`即可执行，并在 odpscmd 中使用 `show tables;` 列出项目中已经存在的表，用 `desc na12878_fastq_merge;` 查看 `na12878_fastq_merge` 这张表的 schema，用 `read na12878_fastq_merge 10;` 来查看这张表的前 10 行，用`CREATE TABLE test_table AS SELECT * FROM fq_r0140900861 limit 10000;` 查看一个表（如`fq_whyd15030427_a`）的记录数，用`select count(*) from fq_whyd15030427_a;`，从表（如fq_r0140900861）中抽取前10000行生成至一个新表`test_table`。
>更多 odpscmd 的使用方法，可以参考我们的官方文档 https://help.aliyun.com/document_detail/27971.html

因为 fastq 的特殊性，使用 odpscmd 中的通用上传命令无法做到 schema 对应，因此在上传数据的时候使用了专门用于做 ETL 的工具 [tattoo](http://repo.aliyun.com/etl)，位于流程中的 `tattoo` 目录。

> **以上，odpscmd 和 tattoo 这两个工具依赖 java 环境，请事先在运行环境中准备好 java >= 1.7**。

MaxCompute账号配置信息 `conf/odps_config.ini` 和 `odpscmd/conf/odps_config.ini` 中（这两个文件相同）。

为了更加方便使用，目前已经将很多常用操作直接封装在流程主程序 `odps_gene.pl` 里，这个脚本会根据参数不同自动调用相关工具完成诸如建表，上传数据，跑不同阶段任务等操作。可以通过阅读 `odps_gene.pl` 以及 `bin/`

### 集群环境准备

- 上传静态资源文件，如: `human_g1k_v37.tar.gz`, `dbsnp_138.tar.gz`, `1000G_phase1.indels.tar.gz`, `Mills_and_1000G_gold_standard.indels.tar.gz`。
注意每个压缩包里应包括基因参考数据和对应的索引文件，具体如下所示。

```
1000G_phase1.indels.tar.gz
-rw-rw-r-- admin/admin 237684308 2016-04-25 10:56:44 1000G_phase1.indels.b37.vcf
-rw-rw-r-- admin/admin   1238166 2016-04-28 17:35:41 1000G_phase1.indels.b37.vcf.idx

dbsnp_138.tar.gz
-rw-rw-r-- admin/admin 10641623137 2016-04-25 11:00:13 dbsnp_138.b37.vcf
-rw-rw-r-- admin/admin    12380770 2016-05-03 16:59:58 dbsnp_138.b37.vcf.idx

human_g1k_v37.tar.gz
-rw-rw-r-- admin/admin   11097 2016-04-22 10:23:43 human_g1k_v37.dict
-rw-rw-r-- admin/admin 3153506519 2016-04-21 18:14:57 human_g1k_v37.fasta
-rw-rw-r-- admin/admin       6597 2016-04-22 09:24:56 human_g1k_v37.fasta.amb
-rw-rw-r-- admin/admin       6844 2016-04-22 09:24:55 human_g1k_v37.fasta.ann
-rw-rw-r-- admin/admin 3101804844 2016-04-22 09:24:00 human_g1k_v37.fasta.bwt
-rw-rw-r-- admin/admin       2746 2016-04-22 10:23:23 human_g1k_v37.fasta.fai
-rw-rw-r-- admin/admin  775451186 2016-04-22 09:24:51 human_g1k_v37.fasta.pac
-rw-rw-r-- admin/admin 1550902424 2016-04-22 09:51:13 human_g1k_v37.fasta.sa

Mills_and_1000G_gold_standard.indels.tar.gz
-rw-rw-r-- admin/admin 86369975 2016-04-22 18:05:13 Mills_and_1000G_gold_standard.indels.b37.vcf
-rw-rw-r-- admin/admin  2356127 2016-04-28 17:19:07 Mills_and_1000G_gold_standard.indels.b37.vcf.idx
```

把这些压缩包放到一个目录下，上传命令如下:

```
./odps_gene.pl -setup ${DIRECTORY}
```

### 数据准备（NA12878为例）

假设输入数据由两个 fastq 文件组成，分别是 `NA12878_1.fastq` 和 `NA_12878_2.fastq`，准备上传到 `bgi_na12878_fastq` 这张 ODPS 表里，命令如下:

```
./odps_gene.pl -upload-fastq NA12878_1.fastq NA12878_2.fastq bgi_na12878_fastq conf/odps_config.ini 8
```

如果某个人的输入文件有多个，可以用不同输入文件对（file-pair），相同输出表重复上述命令，将这个人的多个文件都集中到同一张表里。

## 作业流程

如附件 `全基因组测序.pptx` 所示，数据上传完毕之后，在 ODPS 的执行流程共分 4 个步骤，具体如下：

### 计算每个样本的 gvcf

完成 fastq QC 并调用 bwa、samtools 以及 gatk 完成 Haplotype 的步骤。

> ODPS 中分布式运行的脚本: `mapper_custompart.pl`, `reducer_gvcf.pl`

命令: `./odps_gene.pl -haplotype 输出表 输入表 reducer个数 ${GROUP_HEADER}`

注意事项:
- 参数顺序中**输出**表在前，**输入**表在后
- 899 这个数字应和 bin/part.tsv 的分片数量（行数 `wc -l bin/part.tsv`）一致
- GROUP_HEADER 这个参数从命令行输入时建议用单引号，防止 shell 转义的影响，如 ''
- 此作业的运行时间较长。作业运行时屏幕上会出现形如这样的链接 http://logview.odps.aliyun.com/logview/?h=http://service.odps.aliyun.com/api&p=bgi_speedwgs&i=20161230035845154gzl5qqet&token=TkNzMExRcFNRUHZaL2dQNlRENWh2enlQTXFBPSxPRFBTX09CTzoxNzQzOTM5Nzg4MTM0MDU3LDE0ODM2NzUxMjcseyJTdGF0ZW1lbnQiOlt7IkFjdGlvbiI6WyJvZHBzOlJlYWQiXSwiRWZmZWN0IjoiQWxsb3ciLCJSZXNvdXJjZSI6WyJhY3M6b2RwczoqOnByb2plY3RzL2JnaV9zcGVlZHdncy9pbnN0YW5jZXMvMjAxNjEyMzAwMzU4NDUxNTRnemw1cXFldCJdfV0sIlZlcnNpb24iOiIxIn0= 这是 ODPS 作业的 logview，可以在浏览器中打开这个 url 来随时监控 ODPS 作业的运行情况。作业一旦已经提交（logview 已经生成），即使 **执行命令的终端** 崩溃，进程被杀等，**都不会影响** ODPS 服务端作业的运行。

>关于 logview 的更多信息可以参考官方文档 https://help.aliyun.com/document_detail/27987.html

输入表为上传的 `na12878_fastq`，输出表 `na12878_gvcf`，命令行样例：
```
./odps_gene.pl -haplotype na12878_gvcf na12878_fastq 899 '@RG\tID:HiseqEAAAGAAA-98\tPL:illumina\tPU:150430_I00137_FCH2JTWBBXX_L4_HiseqEAAAGAAA-98\tLB:HiseqEAAAGAAA-98\tSM:NA12878\tCN:GIAB'
```

### 合并所有 gvcf 到一张表

命令: `./odps_gene.pl -union 输出表 输入表1,输入表2,输入表3...`

注意事项:
- 参数顺序中输出表在前，输入表在后
- 输入表之间用半角逗号分割，但是不能有空格
- union 在一起的输入表需要使用相同的分片文件进行计算，如果某些输入表是使用不同的分片文件进行计算，这一步不会报错，但是可能导致后续的计算结果不准或直接出错

假设 na12878，na12891，na12892，na24631 已经完成了上一步的 gvcf 计算，将这些输出表合并的的命令行样例:
```
./odps_gene.pl -union bgi_gvcf_merged na12878_gvcf,na12891_gvcf,na12892_gvcf,na24631_gvcf
```

### 群体 call variant 得到 vcf 表（GenotypeGVCFs）

按照 part.tsv 中指定的分片，对分片内所有样本的 gvcf 调用 gatk 的 genotype 功能进行群体 call 变异，得到分片的 vcf 数据。

ODPS 中分布式运行的脚本: `mapper_genotype.pl`, `reducer_genotype.pl`

命令: `./odps_gene.pl -genotype 输出表 输入表 899`

注意事项:
- 参数顺序中输出表在前，输入表在后
- 899 这个数字应和 bin/part.tsv 的分片数量（行数 `wc -l bin/part.tsv`）一致

样例:
```
./odps_gene.pl -genotype bgi_vcf_for_combine bgi_gvcf_merged 899
```

### 合并 vcf（CombineVariants）

用 gatk 的 CombineVariants 功能把各分片的 vcf 合并成一个完整的 vcf。

ODPS 中分布式运行的脚本: `mapper_combinevcf.pl`, `reducer_combinevcf.pl`

命令: `./odps_gene.pl -combine-vcf 输出表 输入表`

注意事项:
- 只有经过本步骤计算的输出表才能使用 `bin/download.pl` 下载成 vcf 文件，之前的步骤的 gvcf 表等无法被正确下载。

样例:
```
./odps_gene.pl -combine-vcf bgi_vcf bgi_vcf_for_combine
```

## 下载 vcf

脚本 `bin/download.pl` 下载最终的 vcf 表数据到本地，并多分片的排序及还原 header，最终产出标准的 vcf 格式文件。

命令: `bin/download.pl 输入表 输出文件 ODPS配置文件`

注意事项:
- 参数中输入表在前，输出文件在后（跟 `odps_gene.pl` 的风格不同）
- 根据经验，下载的 vcf 文件约 1G，下载以及排序会花一些时间
- 下载以及后续的排序过程可能在当前目录产生一些临时文件，是正常现象，需要保证本地路径有足够硬盘空间。临时文件在流程正常结束后会被自动清除

样例:
```
bin/download.pl bgi_vcf bgi.vcf conf/odps_config.ini
```

## 运算结果简单对比

下载并安装 [vcftools](https://vcftools.github.io/)，对比 vcf 一致性的命令如下:

```
vcftools --vcf result1.vcf --diff result2.vcf --diff-site
```


## 具体执行实例

以下例子是基于已经完整可访问的数据准备的：

### step1 fastq->gvcf

```
time perl odps_gene.pl -haplotype gvcf_r0140900815_3 fq_r0140900815 899 "@RG\tID:HiseqEAAAGAAA-98\tPL:illumina\tPU:150430_I00137_F\tLB:HiseqAAGAAA-98\tSM:R0140900815\tCN:GIAB" > step1.R0140900815.o.log 2> step1.R0140900815.e.log && echo " ** gvcf_r0140900815_3 done **"
```


### Step2 gvcf->vcf

```
time perl odps_gene.pl -union gvcf_merged_r0140900815 gvcf_r0140900815_3 > m.o.log 2> m.e.log && echo "** gvcf_merged_r0140900815 done **"

time perl odps_gene.pl -genotype vcf_for_combine_r0140900815 gvcf_merged_r0140900815 899 > g.o.log 2> g.e.log && echo "** vcf_for_combine_r0140900815 done **"

time perl odps_gene.pl -combine-vcf r0140900815_vcf vcf_for_combine_r0140900815 > c.o.log 2> c.e.log && echo "** r0140900815.vcf done **"
```

### Download The final VCF file

```
bin/download.pl r0140900815_vcf r0140900815.vcf conf/odps_config.ini
```

















