#time ./odps_gene.pl -haplotype tmp_whyd15030427_a_gvcf fq_whyd15030427_a 899 '@RG\tID:HiseqEAAAGAAA-98\tPL:illumina\tPU:150430_I00137_F\tLB:HiseqAAGAAA-98\tSM:WHYD15030427-A\tCN:BGI'  # real-time      142m41.780s

#time ./odps_gene.pl -union tmp_whyd15030427_a_gvcf_comb tmp_whyd15030427_a_gvcf # 1m5.582s
#time ./odps_gene.pl -genotype tmp_whyd15030427_a_genotype_set tmp_whyd15030427_a_gvcf_comb 899  # 16m1.574s
time ./odps_gene.pl -combine-vcf bgi_vcf tmp_whyd15030427_a_genotype_set   # 44m20.009s

