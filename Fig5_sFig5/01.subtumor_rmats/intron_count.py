import pysam
import pandas as pd
import collections
import sys
from multiprocessing import Pool

def get_read_intron(read):
    cigartuples = read.cigartuples
    base_position = read.pos
    for op, nt in cigartuples:
        if op in match_or_deletion:
            base_position += nt
        elif op == BAM_CREF_SKIP:
            junc_start = base_position
            base_position += nt
            return junc_start,base_position

def get_strand_sign(strand):
    if strand == '+':
        return True
    else:
        return False

def find_intron_from_bam(bam_file,chrom,intron_start,intron_end,strand):
    strand_sig = get_strand_sign(strand)
    res = collections.Counter()
    input_bam = pysam.AlignmentFile(bam_file, "rb")
    for read in input_bam.fetch(chrom,intron_start - 6,intron_end + 5, multiple_iterators=True):
        if 'N' in read.cigarstring:
            read_intron_start,read_intron_end = get_read_intron(read)
            if strand_sig:
                if read_intron_start == intron_start and read_intron_end == intron_end and not read.is_reverse:
                    x = str(read.get_tag('Cx'))
                    y = str(read.get_tag('Cy'))
                    barcode = x + '_' + y
                    res[barcode] += 1
            else:
                if read_intron_start == intron_start and read_intron_end == intron_end and read.is_reverse:
                    x = str(read.get_tag('Cx'))
                    y = str(read.get_tag('Cy'))
                    barcode = x + '_' + y
                    res[barcode] += 1

    intron_df = pd.DataFrame(res.items(),columns = ['barcode','count'])
    intron_name = chrom + ":" + str(intron_start) + ":" + str(intron_end) + strand
    intron_df['intron'] = intron_name
    return intron_df

def check_chrom(bam_file):
    input_bam = pysam.AlignmentFile(bam_file, "rb")
    for i in input_bam.references:
        if 'chr' in i:
            return True
    input_bam.close()
    return False

def process_sj_row(intron):
    chrom, intron_start, intron_end = intron.split(':')
    if not chrom_with_chr:
        chrom = chrom.split('chr')[1]
    intron_start = int(intron_start)
    strand = intron_end[-1]
    intron_end = intron_end[:-1]
    intron_end = int(intron_end)
    
    sj_df = find_intron_from_bam(bam_file, chrom, intron_start, intron_end,strand)
    return sj_df

bam_file = sys.argv[1]
intron_file = sys.argv[2]
out_intron_df = sys.argv[3]
intron_ri = sys.argv[4]
out_ri_count = sys.argv[5]
threads = int(sys.argv[6])

match_or_deletion = {0, 2, 7, 8}
BAM_CREF_SKIP = 3
chrom_with_chr = check_chrom(bam_file)
intron_df = pd.read_csv(intron_file,header =None,names = ['intron'])
intron_ri = pd.read_csv(intron_ri,header =None,names = ['intron'])

intron_list = intron_df['intron'].tolist()
intron_ri_list = intron_ri['intron'].tolist()

pool = Pool(processes = threads)
results = []
for sj_df in pool.imap(process_sj_row, intron_list, chunksize = 100):
    results.append(sj_df)
pool.close()
pool.join()

intron_df = pd.concat(results,axis = 0)
intron_df.to_csv(out_intron_df,index=None,compression='gzip')


def process_ri_row(intron):
    chrom, intron_start, intron_end = intron.split(':')
    if not chrom_with_chr:
        chrom = chrom.split('chr')[1]
    intron_start = int(intron_start)
    strand = intron_end[-1]
    intron_end = intron_end[:-1]
    intron_end = int(intron_end)
    
    sj_df = find_ri_from_bam(bam_file, chrom, intron_start, intron_end,strand)
    return sj_df

def find_ri_from_bam(bam_file,chrom,intron_start,intron_end,strand):
    strand_sig = get_strand_sign(strand)
    res = collections.Counter()
    input_bam = pysam.AlignmentFile(bam_file, "rb")
    for read in input_bam.fetch(chrom,intron_start - 1,intron_end, multiple_iterators=True):
        if 'N' in read.cigarstring:
            read_intron_start,read_intron_end = get_read_intron(read)
            if strand_sig:
                if not read.is_reverse:
                    if read_intron_start <= intron_start and intron_end <= read_intron_end:
                        pass
                    else:
                        x = str(read.get_tag('Cx'))
                        y = str(read.get_tag('Cy'))
                        barcode = x + '_' + y
                        res[barcode] += 1
            else:
                if read.is_reverse:
                    if read_intron_start <= intron_start and intron_end <= read_intron_end:
                        pass
                    else:
                        x = str(read.get_tag('Cx'))
                        y = str(read.get_tag('Cy'))
                        barcode = x + '_' + y
                        res[barcode] += 1
        else:
            if strand_sig:
                if not read.is_reverse:
                    x = str(read.get_tag('Cx'))
                    y = str(read.get_tag('Cy'))
                    barcode = x + '_' + y
                    res[barcode] += 1
            else:
                if read.is_reverse:
                    x = str(read.get_tag('Cx'))
                    y = str(read.get_tag('Cy'))
                    barcode = x + '_' + y
                    res[barcode] += 1

    intron_df = pd.DataFrame(res.items(),columns = ['barcode','count'])
    intron_name = chrom + ":" + str(intron_start) + ":" + str(intron_end) + strand
    intron_df['intron'] = intron_name
    return intron_df

pool = Pool(processes = threads)
results = []
for sj_df in pool.imap(process_ri_row, intron_ri_list, chunksize = 100):
    results.append(sj_df)
pool.close()
pool.join()

intron_ri = pd.concat(results,axis = 0)
intron_ri.to_csv(out_ri_count,index=None,compression='gzip')