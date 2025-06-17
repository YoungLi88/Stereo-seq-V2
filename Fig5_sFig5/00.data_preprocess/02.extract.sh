awk '$3 == "exon" {print $1 "\t" $4-1 "\t" $5}' /storage/liuyi/08.stereo_v2/reference/star/human/gencode.v31.chr_patch_hapl_scaff.annotation.gtf > exons.bed
sort -k1,1 -k2,2n exons.bed | bedtools merge -i - > exons_merged.bed
samtools view -b -L exons_merged.bed C03637A5.dup.tumor.bam > C03637A5.dup.tumor.f.bam
samtools view -b -L exons_merged.bed C03637A5.dup.2.bam > C03637A5.dup.2.f.bam
samtools view -b -L exons_merged.bed C03637A5.dup.1.bam > C03637A5.dup.1.f.bam
samtools view -b -L exons_merged.bed C03637A5.dup.normal.bam > C03637A5.dup.normal.bam
samtools index C03637A5.dup.2.f.bam
samtools index C03637A5.dup.1.f.bam
samtools index C03637A5.dup.normal.f.bam
samtools index C03637A5.dup.tumor.f.bam

