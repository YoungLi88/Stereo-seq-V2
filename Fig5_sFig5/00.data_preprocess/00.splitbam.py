import pysam
import os
import gzip
from multiprocess import Pool

def check_bam_index(bam): 
    bai_file = f"{bam}.bai"
    if not os.path.exists(bai_file):
        raise Exception(f"{bam} without index")

def find_chromosomes_from_bam(bam, only_autosome):
    inputfile = pysam.AlignmentFile(bam,"rb")
    chrom_tuple = inputfile.references
    inputfile.close()
    if only_autosome:
        chrom_tuple_tmp = [chrom for chrom in chrom_tuple if chrom.startswith('chr')]
        if len(chrom_tuple_tmp) == 0:
            chrom_tuple_tmp = [chrom for chrom in chrom_tuple if "K" not in chrom and "G" not in chrom and "L" not in chrom and "J" not in chrom and "V" not in chrom and "Q" not in chrom and "N" not in chrom and "h" not in chrom]
        return chrom_tuple_tmp
    else:
        return chrom_tuple
        
def find_chromosomes_from_fasta(fasta_file,only_autosome):
    chromosomes = set()
    if fasta_file.endswith(".gz"):
        f = gzip.open(fasta_file, "rt")
    else:
        f = open(fasta_file)
    for line in f:
        if line.startswith(">"):
            chromosome = line.split()[0][1:]
            chromosomes.add(chromosome)
    f.close()
    return chromosomes
    

def filter_chrom(chromosomes):
    if only_autosome:
        chrom_tuple_tmp = [chrom for chrom in chromosomes if chrom.startswith('chr')]
    if len(chrom_tuple_tmp) == 0:
        chrom_tuple_tmp = [chrom for chrom in chromosomes if "K" not in chrom and "G" not in chrom and "L" not in chrom and "J" not in chrom and "V" not in chrom and "Q" not in chrom and "N" not in chrom and "h" not in chrom]
    po = ';'.join(chrom_tuple_tmp)
    print(f'Find chrom: {po}')
    return chrom_tuple_tmp


def Code_barcode(barcode):
    code = {1:'AA',2:'AT',3:'AC',4:'AG',5:'TA',6:'TT',7:'TC',8:'TG',9:'CA',0:'CC'}
    barcode_list = list()
    for m in barcode:
        barcode_list.append(code[int(m)])
    return ''.join(barcode_list)


def split_bam(bam,chrom,sample_name,output_dir,cluster_map):
    inputfile = pysam.AlignmentFile(bam, "rb")
    
    cluster_list = list(set(cluster_map.values()))
    
    output_dict = {}
    for cluster in cluster_list:
        tmp_bam = os.path.join(output_dir, '%s.%s.%s.bam'%(sample_name,chrom,cluster))
        output = pysam.AlignmentFile(tmp_bam, 'wb', template = inputfile)
        output_dict[cluster] = output
        
    for read in inputfile.fetch(str(chrom), multiple_iterators=True): 
        read_name = read.query_name
        read_seq = read.query_sequence
        read_flag = read.flag
        mapping_q = read.mapping_quality
        
        CB = read_name.split('|||')[1].split('CB:Z:')[1]
        x = Code_barcode(CB.split('_')[0])
        y = Code_barcode(CB.split('_')[1])
        code_CB = x + '-' + y
        read.set_tag('CB',code_CB)
        
        x = int(CB.split('_')[0])
        y = int(CB.split('_')[1])
        binx = x//50*50
        biny = y//50*50
        loc = str(binx) + '_' + str(biny)
        cluster = cluster_map.get(loc,False)
        if cluster:
            output = output_dict[cluster]
            output.write(read)

fasta_file = None
bam = '/media/transport/TNBC/C03637A5_align.tag.bam'
sampleid = "C03637A5"
out_dir = "/storage/liuyi/08.stereo_v2/HER2_breast/bam"
threads = 20 


import scanpy as sc
h5ad_file = "./new.cnv.h5ad"
adata = sc.read_h5ad(h5ad_file)
cluster_map = dict(zip(adata.obs['x'].astype(str) + '_' + adata.obs['y'].astype(str),adata.obs['top_region_cluster']))

check_bam_index(bam)
if fasta_file != None:
    chrom_tuple = find_chromosomes_from_fasta(fasta_file,only_autosome = True)
else:
    chrom_tuple = find_chromosomes_from_bam(bam,only_autosome = True)

split_pool = Pool(threads)
for chrom in chrom_tuple:
    split_pool.apply_async(split_bam, args=(bam,chrom,sampleid,out_dir,cluster_map))
split_pool.close()
split_pool.join()
