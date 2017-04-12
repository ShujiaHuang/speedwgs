# 以下加重点的是需要开放的参数
# BWA MEM 参数
bwa mem (options) *idxbase* *in1.fq* *in2.fq*
options:
## -R str  read group header line such as “@RG\tID:HiseqEAAAGAAA-98\tPL:illumina\t” (null)
-M int  mark shorter split hits as secondary
-t int  number of threads (1)

# samtools 参数
samtools sort in.bam
-@  Set number of sorting and compression threads (1)
-m  Set maximum memory per thread; suffix K/M/G recognized (768M)

samtools rmdup input.sort.bam output.bam
-s  rmdup for SE reads
-S  treat PE reads as SE in rmdup (force -s)

=============================================================================

# GATK 3.5 流程以及参数：
```
gatk -nt 4  -T RealignerTargetCreator -R $refgenomefasta \
	[-l OFF] \
    -I $inbam \
    -o sample.intervals

gatk -T IndelRealigner \
    -R reference.fasta \
    [-l OFF] \
    -I input.bam \
    -rf NotPrimaryAlignment \
    -known Mills_and_1000G_gold_standard.indels.b37.vcf\
	-known 1000G_phase1.indels.b37.vcf\
    -targetIntervals sample.intervals \
    -o realigned.bam

gatk -nct 4 -T BaseRecalibrator \
	[-l OFF] \
	--disable_indel_quals \
	-R $refgenomefasta \
    -knownSites dbsnp_138.b37.vcf \
    -knownSites 1000G_phase1.indels.b37.vcf \
    -knownSites Mills_and_1000G_gold_standard.indels.b37.vcf \
    -I realigned.bam \
    -o recalibration_report.grp

gatk –T PrintReads \
	[-l OFF] \
	-R reference.fasta \
	-BQSR recalibration_report.grp \
	-I realiged.bam \
	-o recalibrated.bam
```

# 下面使用HaplotypeCaller+GVCF的形式，有两种做法，两种都应该尝试，看看结果会不会区别，
我这边线下测试是没区别的：

第一种，先对每个样本的相同区域（如 reducer中规定的那些）各自生成区域性的GVCF文件，
然后Call这些同区域的GVCF，得到这个区域的vcf结果,最终再合并：
=======================================================================
# 生成gvcf
```
for id_start_end in `cat reducer.tsv` ; do
	for i in {1..sample_num}; do 
		gatk -T HaplotypeCaller \
			[-l OFF] \
			#-variant_index_type LINEAR 
			#-variant_index_parameter 128000
			-R $refgenomefasta \
			-I recalibrated.bam \
			--emitRefConfidence GVCF \      # 这里规定输出gvcf
			-L id:start-end \               # 这里规定只在特定的interval中计算
			-o sample${i}.${id_start_end}.g.vcf  # `$id_start_end` 名字中记录这个区间信息，方便处理，也可以用1，2，3,这样的数字代替。
	done

	gatk -T GenotypeGVCFs \
		-R $refgenomefasta \
		--variant sample1.${id_start_end}.g.vcf \
		--variant sample2.${id_start_end}.g.vcf \
		--variant sample3.${id_start_end}.g.vcf \
		...
		--variant sample_n.${id_start_end}.g.vcf \
		-o samples.id_start_end.vcf && echo "done"
done
```

最后合并：
```
gatk -T CombineVariants  \
	-R $fasta \
 	--variant samples.id_start_end1.vcf \
 	--variant samples.id_start_end2.vcf  \
 	...
 	--variant samples.id_start_end2.vcf \
 	-o samples.vcf \
```


第二种，对每个样本各自出完整的gvcf，最后再合并这些不同的样本一起call最后的变异VCF
```
for id_start_end in `cat reducer.tsv` ; do
	gatk -T HaplotypeCaller \
		[-l OFF] \
		-R $refgenomefasta \
		-I recalibrated.bam \
		--emitRefConfidence GVCF \      # 这里规定输出gvcf
		-L id:start-end \               # 这里规定只在固定的interval中计算
		-o sample1.${id_start_end}.g.vcf  # `$id_start_end` 
done

gatk -T CombineGVCFs \
	-R $fasta \
 	--variant sample1.id_start_end1.g.vcf  \
 	--variant sample1.id_start_end2.g.vcf  \
 	...
 	--variant sample1.id_start_end_n.g.vcf \
 	-o samples1.g.vcf \

gatk -T GenotypeGVCFs \
		-R $refgenomefasta \
		--variant sample1.g.vcf \
		--variant sample2.g.vcf \
		--variant sample3.g.vcf \
		...
		--variant sample_n.g.vcf \
		-o samples.vcf && echo "done"

```
************

# 这是最后一步，还有些问题，我还在测试，暂时不跑这一步，需要用上所有的VCF数据。

```
gatk -T VariantRecalibrator \
	-R $refgenomefasta \
	-input sample.vcf \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf \
	-resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.vcf \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf \
	-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 dbsnp_138.b37.vcf 
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
	-recalFile sample.recalibrator \
	-tranchesFile sample_snps.tranches \
	-rscriptFile $sample_recalibrator_plots.R
```

