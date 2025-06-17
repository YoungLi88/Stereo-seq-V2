import gzip

#Input and output file paths
gtf_file = "/genome/10x/mm10/refdata-gex-mm10-2020-A/genes/genes.gtf"
bed_file = "exons.bed"


with open(gtf_file, "r") as gtf, open(bed_file, "w") as bed:
    for line in gtf:
        if line.startswith("#"):
            continue
        
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue
        
        if fields[2] != "exon":
            continue
        # Extract chromosome, start position, end position and chain information
        chrom = fields[0]
        start = int(fields[3]) - 1  #Convert to 0-based
        end = int(fields[4])
        strand = fields[6]
        
        attributes = fields[8]
        gene_id = "unknown"
        
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1]
                break
        
        bed.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t.\t{strand}\n")

