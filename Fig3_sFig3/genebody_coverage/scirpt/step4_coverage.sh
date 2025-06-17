geneBody_coverage=/home/zhufan/anaconda3/envs/cellphonedb/bin/geneBody_coverage.py
gtf_bed="gtf_bed.bed"
clean_bam="/data_hub1/zhufan/02.v2/gc/stereo_v1/CRR1252509/CRR1252509_Aligned_sort.out.bam"
sample_name="CRR1252509"

$geneBody_coverage -r $gtf_bed -i $clean_bam -o $sample_name
