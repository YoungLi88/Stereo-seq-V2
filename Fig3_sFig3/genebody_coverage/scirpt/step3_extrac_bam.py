import pysam
import sys
import pandas as pd
sample = sys.argv[1]
#Reading BAM Files
bamfile = pysam.AlignmentFile(f"{sample}.bam", "rb")

# Read the sorted GC ratio file
df = pd.read_csv('sorted_gene_gc_content.txt', sep='\t')
df['gene'] = df['gene'].apply(lambda x : x.split("(")[0])
# Calculating quantiles
top_10 = df.head(int(len(df) * 0.1))
bottom_10 = df.tail(int(len(df) * 0.1))
middle_10 = df.iloc[int(len(df) * 0.45):int(len(df) * 0.55)]

#Save grouping results
top_10.to_csv('top_10_gc_genes.txt', sep='\t', index=False)
middle_10.to_csv('middle_10_gc_genes.txt', sep='\t', index=False)
bottom_10.to_csv('bottom_10_gc_genes.txt', sep='\t', index=False)

top_genes = set(pd.read_csv('top_10_gc_genes.txt', sep='\t')['gene'])
middle_genes = set(pd.read_csv('middle_10_gc_genes.txt', sep='\t')['gene'])
bottom_genes = set(pd.read_csv('bottom_10_gc_genes.txt', sep='\t')['gene'])

top_bam = pysam.AlignmentFile(f"{sample}_top_10_gc.bam", "wb", template=bamfile)
middle_bam = pysam.AlignmentFile(f"{sample}_middle_10_gc.bam", "wb", template=bamfile)
bottom_bam = pysam.AlignmentFile(f"{sample}_bottom_10_gc.bam", "wb", template=bamfile)
#Iterate over BAM files and split
for read in bamfile:
    #Filter conditions: mapping quality > 10 and not duplicate reads
    if read.mapping_quality > 10 and not read.is_duplicate:
        # 从 Extract gene names from GN tags
        gene_name = read.get_tag('GX') if read.has_tag('GX') else None
        #Group by gene name
        if gene_name in top_genes:
            top_bam.write(read)
        elif gene_name in middle_genes:
            middle_bam.write(read)
        elif gene_name in bottom_genes:
            bottom_bam.write(read)

# close file
bamfile.close()
top_bam.close()
middle_bam.close()
bottom_bam.close()
