import pandas as pd

df = pd.read_csv('gene_gc_content.txt', sep='\t')

# Sort by GC ratio 
df = df.sort_values(by='gc_pct', ascending=False)

df.to_csv('sorted_gene_gc_content.txt', sep='\t', index=False)

