#from Bio.SeqUtils import GC"

from Bio.SeqUtils import gc_fraction

input_fasta = "merged_exons.fasta"
output_gc = "gene_gc_content.txt"

with open(input_fasta, "r") as fasta, open(output_gc, "w") as gc_file:
    #Write header
    gc_file.write("gene\tgc_pct\n")
    
    # Reading FASTA files
    current_gene = None
    sequence = ""
    for line in fasta:
        if line.startswith(">"):
            # If it is a new gene, calculate the GC content of the previous gene
            if current_gene:
                gc_pct = gc_fraction(sequence) * 100  # Convert to percentage
                gc_file.write(f"{current_gene}\t{gc_pct:.2f}\n")
            
            # Extract the current gene name
            current_gene = line.strip()[1:]  # Remove ">"
            sequence = ""  #Reset sequence
        else:
            #Splicing sequence
            sequence += line.strip()
    
    #Calculate the GC content of the last gene
    if current_gene:
        gc_pct = gc_fraction(sequence) * 100  #Convert to percentage
        gc_file.write(f"{current_gene}\t{gc_pct:.2f}\n")

