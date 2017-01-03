# 以下加重点的是需要开放的参数
# BWA MEM 参数
bwa mem (options) *idxbase* *in1.fq* *in2.fq*
options:
```
## -R str  read group header line such as “@RG\tID:HiseqEAAAGAAA-98\tPL:illumina\tSM:NA12878\tCN:GIAB” 
-M int  mark shorter split hits as secondary
-t int  number of threads (1)
```

# samtools 参数
```
samtools sort in.bam
-@  Set number of sorting and compression threads (1)
-m  Set maximum memory per thread; suffix K/M/G recognized (768M)

samtools rmdup input.sort.bam output.bam
-s  rmdup for SE reads
-S  treat PE reads as SE in rmdup (force -s)
```


# GATK 流程以及参数：
```
# vf=10G
gatk -nt 4 --fix_misencoded_quality_scores -T RealignerTargetCreator \ 
    [-l OFF] \
    -R human_g1k_v37.fasta \ 
    -I NA12878_rmdup.bam \ 
    -known 1000G_phase1.indels.b37.vcf \ 
    -known Mills_and_1000G_gold_standard.indels.b37.vcf \
    -o ALN.intervals

gatk -T IndelRealigner \ 
    [-l OFF] \
    -R human_g1k_v37.fasta \ 
    -I NA12878_rmdup.bam \ 
    -known 1000G_phase1.indels.b37.vcf \ 
    -known Mills_and_1000G_gold_standard.indels.b37.vcf \ 
    -o NA12878_realign.bam \ 
    --targetIntervals ALN.intervals

gatk -nt 4 -T BaseRecalibrator -R $refgenomefasta \
    [-l OFF] \
    -knownSites 1000G_phase1.indels.b37.vcf \
    -knownSites Mills_and_1000G_gold_standard.indels.b37.vcf \
    -knownSites dbsnp_138.All.fix.vcf \
    -I NA12878_realign.bam \
    -o sorted.realn.bam.recalibration_report.grp

gatk -nt 4 -T PrintReads -R $refgenomefasta \
    [-l OFF] \
    -I NA12878_realign.bam \
    -BQSR sorted.realn.bam.recalibration_report.grp \
    -o NA12878_realign_BQSR4.bam

```

# 变异检测

关于这一步骤中的HaplotypeCaller有两种方法

（1）直接群体进行HaplotypeCaller，这适合于单样本或者不需要扩展和改动的多样本形式；
（2）每个样本先各自生成gvcf(HaplotypeCaller+GVCF)，然后再进行群体genotype，这种方法其实更加灵活，**推荐！**

## 直接群体进行HaplotypeCaller

```
# haplotype caller by chromosome
gatk  -T HaplotypeCaller \ 
    -R human_g1k_v37.fasta \
    -I NA12878_realign_BQSR4.bam \
    -D dbsnp_138.b37.vcf \ 
    -L Y \
    -o NA12878.HC.X.raw.vcf \ 
    -stand_call_conf 50 \ 
    -A QualByDepth \ 
    -A RMSMappingQuality \ 
    -A MappingQualityRankSumTest \ 
    -A ReadPosRankSumTest \ 
    -A FisherStrand \ 
    -A StrandOddsRatio \ 
    -A Coverage && echo “NA12878.HC.X.raw.vcf done”

# combine vcf
java -jar GenomeAnalysisTK.jar \
   -T CombineVariants \
   -R reference.fasta \
   -genotypeMergeOptions UNIQUIFY
   --variant NA12878.HC.1.raw.vcf \
   --variant NA12878.HC.2.raw.vcf \
   …
   --variant NA12878.HC.Y.raw.vcf \
   --variant NA12878.HC.MT.raw.vcf
   -o NA12878.HC.raw.vcf && echo “CombineVariants done”
```

## 下面使用HaplotypeCaller+GVCF的形式，这里也有有两种做法，两种都应该尝试，看看结果会不会区别，
我这边线下测试是没区别的：

### 第一种，先对每个样本的相同区域（如 reducer中规定的那些）各自生成区域性的GVCF文件，

然后分别Call这些同区域的GVCF，得到这个区域的vcf结果,最终再合并：

* 生成gvcf(这是非常正确的)
```
for id_start_end in `cat reducer.tsv` ; do
	for i in {1..sample_num}; do 
		gatk -T HaplotypeCaller \
			[-l OFF] \
			-variant_index_type LINEAR 
			-variant_index_parameter 128000
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

#【问】同个样品如果分成如上多个区域来跑的时候，是否会被当做多个样本？
#【答】不会。在生成gvcf文件的过程中，gatk会根据@RG（也就是bwa比对的read group）信息中SM（也就是sample name）来判断该份gvcf是否来自同一个样本，如果名字相同，那么就会被认为是同一个样本的，不会产生多样本问题

```

* 结果合并
```
gatk -T CombineVariants  \
	-R $fasta \
 	--variant samples.id_start_end1.vcf \
 	--variant samples.id_start_end2.vcf  \
 	...
 	--variant samples.id_start_end2.vcf \
 	-o samples.vcf \
```


#### 第二种，对每个样本各自出完整的gvcf，最后再合并这些不同的样本一起call最后的变异VCF （慢）
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

## 变异质控

这是最后一步，也是通用的一步，最后得到候选变异结果（vcf format）之后都需要进行的质控，但考虑到速度，MaxCompute中不需要进行。

```
## SNP Recalibrator
java -Xmx4g -jar GenomeAnalysisTK.jar \
   -T VariantRecalibrator \
   -R reference.fasta \
   -input NA12878.HC.raw.vcf \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.vcf \ 
   -resource:omini,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.vcf \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.b37.vcf \ 
   -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 dbsnp_138.b37.vcf \ 
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \ 
   -mode SNP \ 
   -recalFile NA12878.HC.snps.recal \
   -tranchesFile NA12878.HC.snps.tranches \ 
    -rscriptFile NA12878.HC.snps.plots.R \ && echo "VariantRecalibrator done”

java -Xmx4g -jar GenomeAnalysisTK.jar -T ApplyRecalibration \
   -R  human_g1k_v37.fasta \
   -input NA12878.HC.raw.vcf \ 
   --ts_filter_level 99.5 \ 
   -tranchesFile NA12878.HC.snps.tranches \ 
   -recalFile NA12878.HC.snps.recal \
   -mode SNP \
   -o NA12878.HC.snps.filtered.vcf \ && echo "ApplyRecalibration done"

## Indel Recalibrator
java -Xmx4g -jar GenomeAnalysisTK.jar -T VariantRecalibrator \
   -R  human_g1k_v37.fasta \
   -input NA12878.HC.snps.filtered.vcf \
   -resource:mills,known=true,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf \
   -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
   -mode INDEL \
   -recalFile NA12878.HC.snps.indels.recal \
   -tranchesFile NA12878.HC.snps.indels.tranches \
   -rscriptFile NA12878.snps.indels.plots.R \ && echo "VariantRecalibrator done”

time java –Xmx4g -jar GenomeAnalysisTK.jar -T ApplyRecalibration \ 
   -R human_g1k_v37.fasta\
   -input NA12878.HC.snps.filtered.vcf \
   --ts_filter_level 99.0 \
   -tranchesFile NA12878.HC.snps.indels.tranches \
   -recalFile NA12878.snps.indels.recal \
   -mode INDEL \
   -o NA12878.HC.snps.indels.filtered.vcf\ && echo "ApplyRecalibration done"
```

