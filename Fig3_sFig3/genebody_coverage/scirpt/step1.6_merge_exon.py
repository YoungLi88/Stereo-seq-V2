from collections import defaultdict

#Input and output file paths
input_fasta = "exons.fasta"
output_fasta = "merged_exons.fasta"

# Store exon sequences for each gene
gene_sequences = defaultdict(str)

# Reading FASTA files
with open(input_fasta, "r") as fasta:
    current_gene = None
    for line in fasta:
        if line.startswith(">"):
            #Extract gene name (extract gene_id from header)
            header = line.strip()
            current_gene = header.split(":")[0][1:]  #Remove ">" and extract gene_id
        else:
            # Add the sequence to the sequence of the corresponding gene
            gene_sequences[current_gene] += line.strip()

#Write the merged FASTA file
with open(output_fasta, "w") as output:
    for gene, sequence in gene_sequences.items():
        output.write(f">{gene}\n{sequence}\n")

