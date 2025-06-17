import spatialAS
import importlib
from multiprocessing import Pool
import pysam
import collections

rmats = '/storage/liuyi/08.stereo_v2/HER2_breast/bam/tumor_normal_dup_rmats/out'
bam_file = "/storage/liuyi/08.stereo_v2/HER2_breast/bam/C03637A5.dup.total.bam"
rmatsresult = spatialAS.pp.read_rmats(rmats,keep_FDR = False)

input_bam = pysam.AlignmentFile(bam_file, "rb")
exon_list = []
for i in ['SE','MXE','A3SS','A5SS']:
    exon = list(rmatsresult[i]['exon'])
    exon_list.extend(exon)
exon_list = list(set(exon_list))

def get_strand_sign(strand):
    if strand == '+':
        return True
    else:
        return False

def process_sj_row(exon):
    chrom = exon.split('@')[0]
    strand = exon.split('@')[1]
    exon_start = int(exon.split('@')[2].split('_')[0])
    exon_end = int(exon.split('@')[2].split('_')[1])
    sj_df = find_exon_from_bam(bam_file, chrom, exon_start, exon_end,strand)
    return sj_df

def find_exon_from_bam(bam_file,chrom,exon_start,exon_end,strand):
    res = collections.Counter()
    exon_pos = list(range(exon_start,exon_end))
    for read in input_bam.fetch(chrom,exon_start - 1,exon_end, multiple_iterators=True):
        strand_sig = get_strand_sign(strand)
        if strand_sig:
            if not read.is_reverse and bool(set(exon_pos) & set(read.positions)):
                barcode = read.query_name.split('|||')[1].split(':')[-1]
                res[barcode] += 1
        else:
            if read.is_reverse and bool(set(exon_pos) & set(read.positions)):
                barcode = read.query_name.split('|||')[1].split(':')[-1]
                res[barcode] += 1
    exon_df = pd.DataFrame(res.items(),columns = ['barcode','count'])
    exon_name = chrom + ":" + str(exon_start) + ":" + str(exon_end) + strand
    exon_df['exon'] = exon_name
    return exon_df

pool = Pool(processes = 20)
results = []
for sj_df in pool.imap(process_sj_row, exon_list, chunksize = 100):
    results.append(sj_df)
pool.close()
pool.join()
exon_total_df = pd.concat(results,axis = 0)
exon_total_df.to_csv('out.exon.tumor_vs_normal.csv',index=None,compression='gzip')

