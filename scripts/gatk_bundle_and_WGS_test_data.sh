#Known datasets: GATK bundle for human b37 reference
# 
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz.md5
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz.md5
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.indels.b37.vcf.gz.md5
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz.md5
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/hapmap_3.3.b37.vcf.gz.md5
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_omni2.5.b37.vcf.gz.md5
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase3_v4_20130502.sites.vcf.gz
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/1000G_phase3_v4_20130502.sites.vcf.gz.md5

# human b37 reference: 使用的时候需要解压并建（bwa）比对的index
wget -c -O human_g1k_v37.fasta.gz ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

## 创建bwa index
bwa index human_g1k_v37.fasta

# Get test fastq data: NA12878
wget -c -O NA12878_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget -c -O NA12878_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz








