import importlib
from multiprocessing import Pool

import pandas as pd
import numpy as np
import re

rmats = '/storage/liuyi/08.stereo_v2/HER2_breast/bam/cluster_dup_rmats/out'
bam_file = "/storage/liuyi/08.stereo_v2/HER2_breast/bam/C03637A5.dup.total.bam"


def tranid2event(tranid):
    chrom = tranid.split(':')[0]
    chrom = chrom.split('chr')[1]
    pattern = ''.join(re.findall(r'[@|]', tranid))
    if pattern == '@@':  # SE
        event = 'SE'
    elif pattern == '@@@':   # MXE
        event = 'MXE'
    elif pattern == '|@': # A5SS
        event = 'A5SS'
    elif pattern == '@|': # A3SS
        event = 'A3SS'
    elif pattern == '@': # RI
        event = 'RI'
    return event

def extract_SE(transid,strand):
    if strand:
        up,mid,down = transid.split(':+@')
    else:
        down,mid,up = transid.split(':-@')
    _,upstreamES,upstreamEE = up.split(':')
    _,downstreamES,downstreamEE = down.split(':')
    _,skexon_start,skexon_end = mid.split(':')
    return [(upstreamEE,skexon_start),(skexon_end,downstreamES),(upstreamEE,downstreamES)]

def extract_MXE(transid,strand):   # ex1  ex1  ex2 ex2
    if strand:
        up,ex1,ex2,down = transid.split(':+@')
    else:
        down,ex2,ex1,up = transid.split(':-@')
        
    _,upstreamES,upstreamEE = up.split(':')
    _,downstreamES,downstreamEE = down.split(':')
    _,E1S,E1E = ex1.split(':')
    _,E2S,E2E = ex2.split(':')
    return [(upstreamEE,E1S),(E1E,downstreamES),(upstreamEE,E2S),(E2E,downstreamES)]
    
        
def extract_A5SS(transid,strand):    # short long
    if strand:
        up,down = transid.split(':+@')
        shortEE,longExonEnd = up.split(':')[2].split('|')
        chrom,flankingES,flankingEE=down.split(':')
        return [(longExonEnd,flankingES),(shortEE,flankingES)]
    else:
        up,down = transid.split(':-@')
        longExonStart,shortES=up.split(':')[2].split('|')
        chrom,flankingES,flankingEE=down.split(':')
        return [(flankingEE,longExonStart),(flankingEE,shortES)]
    
def extract_A3SS(transid,strand):  # short,long
    if strand:
        up,down = transid.split(':+@')
        chrom,flankingES,flankingEE=up.split(':')
        longExonStart=down.split(':')[1].split('|')[0]
        shortES = down.split(':')[1].split('|')[1]
        return [(flankingEE,longExonStart),(flankingEE,shortES)]
    else:
        up,down = transid.split(':-@')
        chrom=up.split(':')[0]
        flankingES=up.split(':')[1]
        shortEE,longExonEnd=down.split(':')[1].split('|')
        return [(longExonEnd,flankingES),(shortEE,flankingES)]
    
    
def extract_RI(transid,strand): 
    if strand:
        up,down = transid.split(':+@')
    else:
        down,up = transid.split(':-@')
    chrom,upES,upEE = up.split(':') 
    chrom,downES,downEE = down.split(':')
    return [(upEE,downES)]

def strand_define(transid,return_sig = True):
    if return_sig:
        if '+' in transid:
            return True
        elif '-' in transid:
            return False
    else:
        if '+' in transid:
            return '+'
        elif '-' in transid:
            return '-'

def tranid2intron_unverified(transid,return_event = False):
    chrom = transid.split(':')[0]
    pattern = ''.join(re.findall(r'[@|]', transid))
    strand = strand_define(transid)
    strand_sig = strand_define(transid,return_sig = False)
    if pattern == '@@':  # SE
        event = 'SE'
        intron_list = extract_SE(transid,strand)
        
    elif pattern == '@@@':   # MXE
        event = 'MXE'
        intron_list = extract_MXE(transid,strand)
        
    elif pattern == '|@': # A5SS
        event = 'A5SS'
        intron_list = extract_A5SS(transid,strand)
    elif pattern == '@|': # A3SS
        event = 'A3SS'
        intron_list = extract_A3SS(transid,strand)
    elif pattern == '@': #RI
        event = 'RI'
        intron_list = extract_RI(transid,strand)
    
    
    intron_name_list = []
    for intron in intron_list:
        intron_name = chrom + ":" + str(int(intron[0])) + ":" + str(int(intron[1])) + strand_sig
        intron_name_list.append(intron_name)
    
    if return_event:
        return event,intron_name_list
    else:
        return intron_name_list
    

def row_has_duplicates(row):
    return len(row) == len(set(row))
    
def se_func(rmats_se):
    def transform_row(x):
        if x["strand"] == "+":
            tran = f"{x['chr']}:{x['upstreamES']}:{x['upstreamEE']}:+@{x['chr']}:{x['exonStart_0base']}:{x['exonEnd']}:+@{x['chr']}:{x['downstreamES']}:{x['downstreamEE']}"
        else:
            tran = f"{x['chr']}:{x['downstreamES']}:{x['downstreamEE']}:-@{x['chr']}:{x['exonStart_0base']}:{x['exonEnd']}:-@{x['chr']}:{x['upstreamES']}:{x['upstreamEE']}"
        return tran
    rmats_se['tran_id'] = rmats_se.apply(transform_row, axis=1)
    return rmats_se

# MXE
def mxe_func(rmats_mxe):
    def transform_row(x):
        if x["strand"] == "+":
            tran = f"{x['chr']}:{x['upstreamES']}:{x['upstreamEE']}:+@{x['chr']}:{x['1stExonStart_0base']}:{x['1stExonEnd']}:+@{x['chr']}:{x['2ndExonStart_0base']}:{x['2ndExonEnd']}:+@{x['chr']}:{x['downstreamES']}:{x['downstreamEE']}"
        else:
            tran = f"{x['chr']}:{x['downstreamES']}:{x['downstreamEE']}:-@{x['chr']}:{x['2ndExonStart_0base']}:{x['2ndExonEnd']}:-@{x['chr']}:{x['1stExonStart_0base']}:{x['1stExonEnd']}:-@{x['chr']}:{x['upstreamES']}:{x['upstreamEE']}"
        return tran
    rmats_mxe['tran_id'] = rmats_mxe.apply(transform_row, axis=1)
    return rmats_mxe

# A5SS
def a5ss_func(rmats_a5ss):
    def transform_row(x):
        if x["strand"] == "+":
            tran = f"{x['chr']}:{x['longExonStart_0base']}:{x['shortEE']}|{x['longExonEnd']}:+@{x['chr']}:{x['flankingES']}:{x['flankingEE']}"
        else:
            tran = f"{x['chr']}:{x['longExonEnd']}:{x['longExonStart_0base']}|{x['shortES']}:-@{x['chr']}:{x['flankingES']}:{x['flankingEE']}"
        return tran

    rmats_a5ss['tran_id'] = rmats_a5ss.apply(transform_row, axis=1)
    return rmats_a5ss

# A3SS
def a3ss_func(rmats_a3ss):
    def transform_row(x):
        if x["strand"] == "+":
            tran = f"{x['chr']}:{x['flankingES']}:{x['flankingEE']}:+@{x['chr']}:{x['longExonStart_0base']}|{x['shortES']}:{x['longExonEnd']}"
        else:
            tran = f"{x['chr']}:{x['flankingES']}:{x['flankingEE']}:-@{x['chr']}:{x['shortEE']}|{x['longExonEnd']}:{x['longExonStart_0base']}"
        return tran

    rmats_a3ss['tran_id'] = rmats_a3ss.apply(transform_row, axis=1)
    return rmats_a3ss

def ri_func(rmats_ri):
    def transform_row(x):
        if x["strand"] == "+":
            tran = f"{x['chr']}:{x['upstreamES']}:{x['upstreamEE']}:+@{x['chr']}:{x['downstreamES']}:{x['downstreamEE']}"
        else:
            tran = f"{x['chr']}:{x['downstreamEE']}:{x['downstreamES']}:-@{x['chr']}:{x['upstreamEE']}:{x['upstreamES']}"
        return tran

    rmats_ri['tran_id'] = rmats_ri.apply(transform_row, axis=1)
    return rmats_ri

# def calculate_position(original_start, original_end, insertion_start, insertion_end, strand = '+'):
#     '''
#     0 → 5‘   1 → 3’
#     '''
#     original_length = abs(original_end - original_start)
#     insertion_mid = (insertion_start + insertion_end) / 2
#     distance_from_start = abs(insertion_mid - original_start)
#     position_percentage = distance_from_start / original_length
#     if strand =='-':
#         position_percentage = 1 - position_percentage
#     return position_percentage

# def compute_se(x):
#     exonStart_0base,exonEnd,strand,gene_start,gene_end = x[['exonStart_0base','exonEnd','strand','gene_start','gene_end']]
#     return calculate_position(gene_start,gene_end,exonStart_0base,exonEnd,strand)

# def compute_mxe(x):
#     es_1,ee_1,es_2,ee_2,strand,gene_start,gene_end = x[['1stExonStart_0base','1stExonEnd','2ndExonStart_0base','2ndExonEnd','strand','gene_start','gene_end']]
#     exon_1 = calculate_position(gene_start,gene_end,es_1,ee_1,strand)
#     exon_2 = calculate_position(gene_start,gene_end,es_2,ee_2,strand)
#     return np.mean([exon_1,exon_2])

# def compute_a3ss(x):
#     les,ses,see,lee,strand,gene_start,gene_end = x[['longExonStart_0base','shortES','shortEE','longExonEnd','strand','gene_start','gene_end']]
#     if strand == '+':
#         return calculate_position(gene_start,gene_end,les,ses,strand)
#     else:
#         return calculate_position(gene_start,gene_end,see,lee,strand)
    
# def compute_a5ss(x):
#     les,ses,see,lee,strand,gene_start,gene_end = x[['longExonStart_0base','shortES','shortEE','longExonEnd','strand','gene_start','gene_end']]
#     if strand == '+':
#         return calculate_position(gene_start,gene_end,see,lee,strand)
#     else:
#         return calculate_position(gene_start,gene_end,les,ses,strand)

def combine_columns(row, sep='_'):
    return sep.join(row.astype(str))

def intron_check(tran_id_df,remain_intron,event,allow_missing = False):
    tran_jug = tran_id_df.isin(remain_intron)
    if allow_missing:
        if event == 'SE':
            allow_tran_id = tran_jug[tran_jug['i3'] & (tran_jug['i1'] | tran_jug['i2'])].index
            tran_id_df = tran_id_df.loc[allow_tran_id].copy()
        elif event == 'MXE':
            allow_tran_id = tran_jug[tran_jug['i2'] | tran_jug['i3']].index
            tran_id_df = tran_id_df.loc[allow_tran_id].copy()
        elif event == 'A3SS' or event == 'A5SS':
            allow_tran_id = tran_jug[tran_jug['i1'] & tran_jug['i2']].index
            tran_id_df = tran_id_df[np.sum(tran_jug,axis = 1) == tran_id_df.shape[1]].copy()
        elif event == 'RI':
            tran_id_df = tran_id_df[np.sum(tran_jug,axis = 1) == tran_id_df.shape[1]].copy()
    else:
        tran_id_df = tran_id_df[np.sum(tran_jug,axis = 1) == tran_id_df.shape[1]].copy()
        
    tran_id_df = tran_id_df.applymap(lambda x:  x  if x in remain_intron else None)
    check_id = tran_id_df.copy()
    check_id['replace_id'] = check_id.apply(combine_columns,axis = 1)
    if event =='A3SS' or event == 'A5SS':
        only_ass = check_id[['i1','i2']]
        check_id = check_id.loc[only_ass[only_ass.apply(row_has_duplicates,axis = 1)].index]
    unique_as_index = check_id.drop_duplicates('replace_id',keep='first').index
    tran_id_df = tran_id_df.loc[unique_as_index].copy()
    return tran_id_df

class as_class:
    def __init__(self,rmats_result):
        self.rmats_result = rmats_result
        self.load_rmats_data()
        self.bulid_intron()
    
    def load_bigwig(self,bigwig):
        self.bigwig = bigwig

    def load_rmats_data(self):
        input_proxy = self.rmats_result
        for event_type in ["SE",'MXE','A3SS','A5SS','RI']:
            file_path = f"{input_proxy}/fromGTF.{event_type}.txt"
            rmats_data = pd.read_csv(file_path, sep='\t')
            if event_type=='SE':
                rmats_data = se_func(rmats_data)
                rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
                self.raw_se = rmats_data                
            elif event_type=='MXE':
                rmats_data = mxe_func(rmats_data)
                rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
                self.raw_mxe = rmats_data
            elif event_type=='A3SS':
                rmats_data = a3ss_func(rmats_data)
                rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
                self.raw_a3ss = rmats_data
            elif event_type=='A5SS':
                rmats_data = a5ss_func(rmats_data)
                rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
                self.raw_a5ss = rmats_data
            elif event_type =='RI':
                rmats_data = ri_func(rmats_data)
                rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
                self.raw_ri = rmats_data


    def select_tranid(self,tranid_list,qcpass=False,scpass=False):
        tmp_se = self.valid_se[self.valid_se['tran_id'].isin(tranid_list)].reset_index(drop=True).copy()
        tmp_mxe = self.valid_mxe[self.valid_mxe['tran_id'].isin(tranid_list)].reset_index(drop=True).copy()
        tmp_a3ss = self.valid_a3ss[self.valid_a3ss['tran_id'].isin(tranid_list)].reset_index(drop=True).copy()
        tmp_a5ss = self.valid_a5ss[self.valid_a5ss['tran_id'].isin(tranid_list)].reset_index(drop=True).copy()
        tmp_ri = self.valid_ri[self.valid_ri['tran_id'].isin(tranid_list)].reset_index(drop=True).copy()
        
        if qcpass:
            print(f'QCpass tran_id:')
            total = 0
            for event,event_df in {'SE':tmp_se,'MXE':tmp_mxe,'A3SS':tmp_a3ss,'A5SS':tmp_a5ss,'RI':tmp_ri}.items():
                print(f'\t{event}:{event_df.shape[0]}')
                total += event_df.shape[0]
            print(f'\tTotal AS event:{total}')
        else:
            total = 0
            for event,event_df in {'SE':tmp_se,'MXE':tmp_mxe,'A3SS':tmp_a3ss,'A5SS':tmp_a5ss,'RI':tmp_ri}.items():
                print(f'\t{event}:{event_df.shape[0]}')
                total += event_df.shape[0]
            print(f'\tTotal AS event:{total}')

        if qcpass:
            self.qcpass_se = tmp_se
            self.qcpass_mxe = tmp_mxe
            self.qcpass_a3ss = tmp_a3ss
            self.qcpass_a5ss = tmp_a5ss
            self.qcpass_ri = tmp_ri
        elif scpass:
            self.scpass_se = tmp_se
            self.scpass_mxe = tmp_mxe
            self.scpass_a3ss = tmp_a3ss
            self.scpass_a5ss = tmp_a5ss
            self.scpass_ri = tmp_ri

    def bulid_intron(self):
        tran_id_se = pd.DataFrame(self.raw_se['tran_id'])
        tran_id_se[['i1','i2','i3']] = tran_id_se['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
        tran_id_se = tran_id_se.set_index('tran_id')

        tran_id_mxe = pd.DataFrame(self.raw_mxe['tran_id'])
        tran_id_mxe[['i1','i2','i3','i4']] = tran_id_mxe['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
        tran_id_mxe = tran_id_mxe.set_index('tran_id')

        tran_id_a3ss = pd.DataFrame(self.raw_a3ss['tran_id'])
        tran_id_a3ss[['i1','i2']] = tran_id_a3ss['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
        tran_id_a3ss = tran_id_a3ss.set_index('tran_id')

        tran_id_a5ss = pd.DataFrame(self.raw_a5ss['tran_id'])
        tran_id_a5ss[['i1','i2']] = tran_id_a5ss['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
        tran_id_a5ss = tran_id_a5ss.set_index('tran_id')

        tran_id_ri = pd.DataFrame(self.raw_ri['tran_id'])
        tran_id_ri[['i']] = tran_id_ri['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
        tran_id_ri = tran_id_ri.set_index('tran_id')
        
        intron_list = pd.concat([
                tran_id_se['i1'],tran_id_se['i2'],tran_id_se['i3'],
                tran_id_mxe['i1'],tran_id_mxe['i2'],tran_id_mxe['i3'],tran_id_mxe['i4'],
                tran_id_a3ss['i1'],tran_id_a3ss['i2'],
                tran_id_a5ss['i1'],tran_id_a5ss['i2'],
                tran_id_ri['i']],axis = 0,ignore_index=True).drop_duplicates()
        intron_in_tranid = pd.DataFrame(intron_list,columns=['intron'])
        
        self.intron_df = intron_in_tranid.reset_index(drop=True)
        self.tmp = [tran_id_se,tran_id_mxe,tran_id_a3ss,tran_id_a5ss,tran_id_ri]

    def check_tranid_intron(self,intron_adata,ri_adata,allow_missing = True):
        
        tran_id_se,tran_id_mxe,tran_id_a3ss,tran_id_a5ss,tran_id_ri = self.tmp
        
        remain_intron = intron_adata.var['gene_ids'].tolist()
        
        tran_id_se = intron_check(tran_id_se,remain_intron,event = 'SE',allow_missing = allow_missing)
        tran_id_mxe = intron_check(tran_id_mxe,remain_intron,event = 'MXE',allow_missing = allow_missing)
        tran_id_a3ss = intron_check(tran_id_a3ss,remain_intron,event = 'A3SS',allow_missing = allow_missing)
        tran_id_a5ss = intron_check(tran_id_a5ss,remain_intron,event = 'A5SS',allow_missing = allow_missing)
        tran_id_ri = intron_check(tran_id_ri,remain_intron,event = 'RI',allow_missing = allow_missing)
        tran_id_ri = tran_id_ri[tran_id_ri['i'].isin(ri_adata.var_names)].copy()
        
        self.tmp = [tran_id_se,tran_id_mxe,tran_id_a3ss,tran_id_a5ss,tran_id_ri]
        
        raw_se = self.raw_se
        raw_mxe = self.raw_mxe
        raw_a3ss = self.raw_a3ss
        raw_a5ss = self.raw_a5ss
        raw_ri = self.raw_ri
        self.valid_se = raw_se[raw_se['tran_id'].isin(tran_id_se.index)].copy()
        self.valid_mxe = raw_mxe[raw_mxe['tran_id'].isin(tran_id_mxe.index)].copy()
        self.valid_a3ss = raw_a3ss[raw_a3ss['tran_id'].isin(tran_id_a3ss.index)].copy()
        self.valid_a5ss = raw_a5ss[raw_a5ss['tran_id'].isin(tran_id_a5ss.index)].copy()
        self.valid_ri = raw_ri[raw_ri['tran_id'].isin(tran_id_ri.index)].copy()
        
        print('>> Write valid tran_id in as_meta.valid_xxxx["se","mxe","a3ss","a5ss"] usage: adata.uns["AS"].tranid["*tran_id"]')
    
        se_sum,mxe_sum,a3ss_sum,a5ss_sum,ri_sum = tran_id_se.shape[0],tran_id_mxe.shape[0],tran_id_a3ss.shape[0],tran_id_a5ss.shape[0],tran_id_ri.shape[0]
        print(f"\tSE: {se_sum}\n\tMXE: {mxe_sum}\n\tA3SS: {a3ss_sum}\n\tA5SS: {a5ss_sum}\n\tRI: {ri_sum}\n\tTotal AS event:{se_sum+mxe_sum+a3ss_sum+a5ss_sum+ri_sum}")
        tranid = {}
        for event,df in {"SE":tran_id_se,"MXE":tran_id_mxe,"A3SS":tran_id_a3ss,"A5SS":tran_id_a5ss,'RI':tran_id_ri}.items():
            tmp_dict = df.apply(lambda row: row.tolist(), axis=1).to_dict()
            tranid.update(tmp_dict)
        self.tranid = tranid
    
    def collect_tranid_meta(self,gtf):
        print('Load GTF From GTF ... It takes times')
        gene_info = pd.read_csv(gtf, sep='\t', comment='#', 
                            usecols=[0,2,3,4,8], names=['chr','type','start','end', 'info'],low_memory=False)
        gene_info = gene_info[gene_info['type'] =='gene']
        gene_info['gene_id'] = gene_info['info'].str.extract(r'gene_id "([^"]+)";')
        start_dict = dict(zip(gene_info['gene_id' ],gene_info['start']))
        end_dict = dict(zip(gene_info['gene_id'],gene_info['end']))
        
        tran_meta = {}
        for event,df in {"SE":self.valid_se,"MXE":self.valid_mxe,"A3SS":self.valid_a3ss,"A5SS":self.valid_a5ss,'RI':self.valid_ri}.items():
            df = df.copy()
            df['gene_start'] = df['GeneID'].map(start_dict)
            df['gene_end'] = df['GeneID'].map(end_dict)
            df['gene_len'] = df['gene_end']- df['gene_start']
            if event == 'SE':
                df['event'] = 'SE'
            elif event == 'MXE':
                df['event'] = 'MXE'
            elif event == 'A3SS':
                df['event'] = 'A3SS'
            elif event == 'A5SS':
                df['event'] = 'A5SS'
            elif event == 'RI':
                df['event'] = 'RI'
            sub_tran2gene_and_pos = {tid: [gl,gs,event] for tid,gl, gs, event in zip(df['tran_id'],df['gene_len'], df['geneSymbol'],df['event'])}
            tran_meta.update(sub_tran2gene_and_pos)
        meta_df = pd.DataFrame(tran_meta).T
        meta_df.columns = ['GeneLenth','GeneSymbol','Event']
        self.tranid_meta = meta_df
        self.tran2genelen = dict(zip(meta_df.index,meta_df['GeneLenth']))
        self.tran2gene = dict(zip(meta_df.index,meta_df['GeneSymbol']))
        self.tran2event = dict(zip(meta_df.index,meta_df['Event']))
        print(f">> Save tranid meta info in AS.tranid_meta")
    
def get_intron_df(rmats_result):
    as_meta = as_class(rmats_result)
    intron_df = as_meta.intron_df
    return intron_df
    
def load_as(adata,ri_adata,rmats_result,gtf,allow_missing = True):
    adata.uns['arg'] = {'allow_missing':allow_missing}
    as_meta = as_class(rmats_result)
    as_meta.check_tranid_intron(adata,ri_adata,allow_missing = allow_missing)
    as_meta.collect_tranid_meta(gtf)
    adata.uns['AS'] = as_meta
    print("Write AS in adata.uns['AS']")
    return adata


def read_rmats(rmats_result,keep_FDR = True):
    rmats_df_dict = {}
    input_proxy = rmats_result
    for event_type in ["SE",'MXE','A3SS','A5SS','RI']:
        file_path = f"{input_proxy}/{event_type}.MATS.JC.txt"
        rmats_data = pd.read_csv(file_path, sep='\t')
        if event_type=='SE':
            rmats_data = se_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            rmats_data['exon'] = rmats_data['chr'] + '@' + rmats_data['strand'] + '@' + rmats_data['exonStart_0base'].map(str) + '_' + rmats_data['exonEnd'].map(str)
        elif event_type=='MXE':
            rmats_data = mxe_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            rmats_data['exon'] = rmats_data['chr'] + '@' + rmats_data['strand'] + '@' + rmats_data['1stExonStart_0base'].map(str) + '_' +  rmats_data['1stExonEnd'].map(str) + '@' + rmats_data['2ndExonStart_0base'].map(str) + '_' + rmats_data['2ndExonEnd'].map(str)
        elif event_type=='A3SS':
            rmats_data = a3ss_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            rmats_data['exon'] = rmats_data['tran_id'].map(lambda x :  x.split(':')[0] + '@' + x.split('@')[0].split(':')[-1] + '@'+ '_'.join(x.split('@')[1].split(':')[1].split('|')))

        elif event_type=='A5SS':
            rmats_data = a5ss_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            rmats_data['exon'] = rmats_data['tran_id'].map(lambda x : x.split(':')[0] + '@' + x.split('@')[0].split(':')[-1] + '@' + '_'.join(x.split('@')[0].split(':')[2].split('|')))

        elif event_type =='RI':
            rmats_data = ri_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            rmats_data['exon'] = rmats_data['chr'] + '@' + rmats_data['strand'] + '@' + rmats_data['upstreamEE'].map(str) + '_' + rmats_data['downstreamES'].map(str)
        if keep_FDR:
            rmats_data = rmats_data[rmats_data['FDR']<=0.05]
        rmats_df_dict[event_type] = rmats_data
    return rmats_df_dict




rmatsresult = read_rmats(rmats,keep_FDR = False)

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
exon_total_df.to_csv('out.exon.tumor_subclone.csv',index=None,compression='gzip')
