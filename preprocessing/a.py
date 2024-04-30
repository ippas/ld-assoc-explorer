import os
import re
import numpy as np
import pandas as pd
import scipy.sparse as sparse

#find ld file by bp chromosome 
def find_ld_file(bp, chrom, folder) :
    for filename in os.listdir(folder):
        if filename.endswith('.gz'):
            match = re.match(r'chr(\d+)_(\d+)_(\d+).gz', filename)
            if match:
                if chrom == int(match.group(1)):
                    begin = int(match.group(2))
                    end = int(match.group(3))
                    if begin <= bp <= end:
                        return filename
    return None

def load_ld_npz(ld_prefix):

    #load the SNPs metadata
    gz_file = '%s.gz'%(ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps['A1'] + '.' + df_ld_snps['A2']

    #load the LD matrix
    npz_file = '%s.npz'%(ld_prefix)
    try:
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s'%(npz_file))

    #create df_R and return it
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps


with open ('../data/snp/chromosome10.csv', 'r') as file:
    for line in file:
        rsid, bp = line.strip().split(',')
        bp = int(bp)
        chr = 10
        ld_filepath = '/slow/projects/ifpan-paula-bbuk'
        ld_file_name = find_ld_file(bp, chr, ld_filepath)

    if ld_file_name:
        print(f'LD file for {rsid} found - {ld_file_name}')
    else:
        print(f'Ld file for {rsid} not found')
