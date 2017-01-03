#time perl odps_gene.pl -union gvcf_merged_r0140900815 gvcf_r0140900815_3 > m.o.log 2> m.e.log && echo "** gvcf_merged_r0140900815 done **"

#time perl odps_gene.pl -genotype vcf_for_combine_r0140900815 gvcf_merged_r0140900815 899 > g.o.log 2> g.e.log && echo "** vcf_for_combine_r0140900815 done **" 

time perl odps_gene.pl -combine-vcf r0140900815_vcf vcf_for_combine_r0140900815 > c.o.log 2> c.e.log && echo "** r0140900815.vcf done **"

