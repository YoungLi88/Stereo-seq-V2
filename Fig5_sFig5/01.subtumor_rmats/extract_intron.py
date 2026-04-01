import sys
rmats_result = sys.argv[1]
outdir = sys.argv[2]
sample = sys.argv[3]

intron_df_file = f'{outdir}/{sample}.Intron.set'
intron_RI_file = f'{outdir}/{sample}.RI.set'

import pandas as pd
import re

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


def load_rmats_data(rmats_result):
    input_proxy = rmats_result
    for event_type in ["SE",'MXE','A3SS','A5SS','RI']:
        file_path = f"{input_proxy}/fromGTF.{event_type}.txt"
        rmats_data = pd.read_csv(file_path, sep='\t')
        if event_type=='SE':
            rmats_data = se_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            raw_se = rmats_data                
        elif event_type=='MXE':
            rmats_data = mxe_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            raw_mxe = rmats_data
        elif event_type=='A3SS':
            rmats_data = a3ss_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            raw_a3ss = rmats_data
        elif event_type=='A5SS':
            rmats_data = a5ss_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            raw_a5ss = rmats_data
        elif event_type =='RI':
            rmats_data = ri_func(rmats_data)
            rmats_data = rmats_data[rmats_data['chr'].isin([f'chr{x}' for x in range(1,23)] + ['chrX','chrY','chrM'])]
            raw_ri = rmats_data
    return raw_se,raw_mxe,raw_a3ss,raw_a5ss,raw_ri

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

def extract_RI(transid,strand):  # short,long
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


raw_se,raw_mxe,raw_a3ss,raw_a5ss,raw_ri = load_rmats_data(rmats_result)
tran_id_se = pd.DataFrame(raw_se['tran_id'])
tran_id_se[['i1','i2','i3']] = tran_id_se['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
tran_id_se = tran_id_se.set_index('tran_id')

tran_id_mxe = pd.DataFrame(raw_mxe['tran_id'])
tran_id_mxe[['i1','i2','i3','i4']] = tran_id_mxe['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
tran_id_mxe = tran_id_mxe.set_index('tran_id')


tran_id_a3ss = pd.DataFrame(raw_a3ss['tran_id'])
tran_id_a3ss[['i1','i2']] = tran_id_a3ss['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
tran_id_a3ss = tran_id_a3ss.set_index('tran_id')


tran_id_a5ss = pd.DataFrame(raw_a5ss['tran_id'])
tran_id_a5ss[['i1','i2']] = tran_id_a5ss['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
tran_id_a5ss = tran_id_a5ss.set_index('tran_id')

tran_id_ri = pd.DataFrame(raw_ri['tran_id'])
tran_id_ri[['i']] = tran_id_ri['tran_id'].apply(lambda x :tranid2intron_unverified(x)).tolist()
tran_id_ri = tran_id_ri.set_index('tran_id')

intron_list = pd.concat([
        tran_id_se['i1'],tran_id_se['i2'],tran_id_se['i3'],
        tran_id_mxe['i1'],tran_id_mxe['i2'],tran_id_mxe['i3'],tran_id_mxe['i4'],
        tran_id_a3ss['i1'],tran_id_a3ss['i2'],
        tran_id_a5ss['i1'],tran_id_a5ss['i2'],
        tran_id_ri['i']
],axis = 0,ignore_index=True).drop_duplicates()
intron_in_tranid = pd.DataFrame(intron_list,columns=['intron'])
intron_df = intron_in_tranid.reset_index(drop=True)

intron_df.to_csv(intron_df_file,header=None,index=None)
tran_id_ri['i'].reset_index()['i'].to_csv(intron_RI_file,header=None,index=None)