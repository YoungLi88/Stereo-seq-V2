#Extract exon sequences
bedtools getfasta   -fi "/genome/10x/mm10/refdata-gex-mm10-2020-A/fasta/genome.fa" -bed exons_sorted.bed -s  -name  -fo exons.fasta
