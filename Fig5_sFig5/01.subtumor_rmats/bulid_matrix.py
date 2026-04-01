import pandas as pd
import numpy as np
import gzip
import os

def InfoFromGem_n(tissue_gem):
    if '.gz' in tissue_gem:
        with gzip.open(tissue_gem, 'rb') as file:
            for line in file:
                line = line.decode('utf-8').strip()
                if line.startswith('#OffsetX='):
                    offset_x = int(line.split('=')[1])
                elif line.startswith('#OffsetY='):
                    offset_y = int(line.split('=')[1])
                    break
    else:
        with open(tissue_gem, 'r') as file:
            for line in file:
                line = line.decode('utf-8').strip()
                if line.startswith('#OffsetX='):
                    offset_x = int(line.split('=')[1])
                elif line.startswith('#OffsetY='):
                    offset_y = int(line.split('=')[1])
                    break    
    return offset_x,offset_y

def InfoFromGem(tissue_gem):
    if '.gz' in tissue_gem:
        with gzip.open(tissue_gem, 'rb') as file:
            for line in file:
                line = line.decode('utf-8').strip()
                if line.startswith('#OffsetX='):
                    offset_x = int(line.split('=')[1])
                elif line.startswith('#OffsetY='):
                    offset_y = int(line.split('=')[1])
                    break
    else:
        with open(tissue_gem, 'r') as file:
            for line in file:
                line = line.decode('utf-8').strip()
                if line.startswith('#OffsetX='):
                    offset_x = int(line.split('=')[1])
                elif line.startswith('#OffsetY='):
                    offset_y = int(line.split('=')[1])
                    break
                    
    gem = pd.read_csv(tissue_gem,sep = '\t',comment='#')
    gem_xmin = gem['x'].min()
    gem_ymin = gem['y'].min()
    return gem_xmin,gem_ymin,offset_x,offset_y

def merge_bin_coor(coor: np.ndarray, coor_min: int, bin_size: int):
    return np.floor((coor - coor_min) / bin_size).astype(np.int32)

def get_bin_center(bin_coor: np.ndarray, coor_min: int, bin_size: int):
    return bin_coor * bin_size + coor_min + int(bin_size / 2)



def tomatrix(outdir_mtx,df):

    if not os.path.exists(outdir_mtx):
        os.makedirs(outdir_mtx)

    df = df.sort_values('barcode')
    barcodes = df['barcode'].drop_duplicates(keep='first').reset_index()
    barcodes_dict = dict(zip(barcodes['barcode'], barcodes.index + 1))
    features = df['feature'].drop_duplicates(keep='first').reset_index()
    features_dict = dict(zip(features['feature'], features.index + 1))
    df['barcode'] = df['barcode'].apply(lambda x: barcodes_dict[x])
    df['feature'] = df['feature'].apply(lambda x: features_dict[x])
    hd = pd.DataFrame([['%%MatrixMarket matrix coordinate integer general'],\
                                   ['%'], \
                                   [' '.join([str(len(features)), str(len(barcodes)), str(len(df))])]])

    hd.to_csv(outdir_mtx+'/matrix.mtx.gz', compression='gzip', sep = '\t', index=False, header=False)
    df[['feature', 'barcode', 'count']].to_csv(outdir_mtx+'/matrix.mtx.gz', compression='gzip', mode='a+',\
                                                                                                     sep=' ', header=False, index=False)
    barcodes['barcode'].to_csv(outdir_mtx+'/barcodes.tsv.gz', compression='gzip', sep='\t', header=False, index=False)
    features['gene'] = features['feature']
    features['feature_type'] = 'intron'
    features[['feature', 'gene', 'feature_type']].to_csv(outdir_mtx+'/features.tsv.gz', compression='gzip', sep='\t', header=False, index=False)

import sys
intron_count=sys.argv[1]
outdir = sys.argv[2]
bin_size = int(sys.argv[3])
gem = sys.argv[4]
bin_type = sys.argv[5]

intron_count = pd.read_csv(intron_count,compression='gzip')
if bin_type == 'normal':
    offset_x,offset_y = InfoFromGem_n(gem)
    split_values = intron_count['barcode'].str.split('_', expand=True)
    intron_count[['bam_x','bam_y']] = split_values[[0,1]].astype(int)
    intron_count['gem_x'] = intron_count['bam_x'] - offset_x
    intron_count['gem_y'] = intron_count['bam_y'] - offset_y
    intron_count['bin_x'] = intron_count['gem_x'] //bin_size * bin_size
    intron_count['bin_y'] = intron_count['gem_y'] //bin_size * bin_size
    intron_count['barcode'] = 'DNB_' + intron_count['bin_x'].map(str) + '_' + intron_count['bin_y'].map(str)
elif bin_type == 'stereo':
    gem_xmin,gem_ymin,offset_x,offset_y = InfoFromGem(gem)
    split_values = intron_count['barcode'].str.split('_', expand=True)
    intron_count[['bam_x','bam_y']] = split_values[[0,1]].astype(int)
    intron_count['gem_x'] = intron_count['bam_x'] - offset_x
    intron_count['gem_y'] = intron_count['bam_y'] - offset_y
    intron_count['bin_x'] = merge_bin_coor(intron_count['gem_x'].values, gem_xmin, bin_size)
    intron_count['bin_y'] = merge_bin_coor(intron_count['gem_y'].values, gem_ymin, bin_size)
    intron_count['x_center'] = get_bin_center(intron_count['bin_x'], gem_xmin, bin_size)
    intron_count['y_center'] = get_bin_center(intron_count['bin_y'], gem_ymin, bin_size)
    intron_count['barcode'] = intron_count['x_center'].map(str) + '_' + intron_count['y_center'].map(str)

ret = intron_count.pivot_table(index=['barcode','intron'], values='count', aggfunc='sum')
ret = ret.reset_index()
ret = ret.rename({'intron':'feature'},axis = 1)
tomatrix(outdir,ret)
