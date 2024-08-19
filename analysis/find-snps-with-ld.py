import os
import re
import sys
import numpy as np
import pandas as pd
from scipy.sparse import load_npz

#function from README in dataset
def load_ld_npz(ld_prefix):

    #load the SNPs metadata
    gz_file = '%s.gz'%(ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2'}, inplace=True)
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps['A1'] + '.' + df_ld_snps['A2']

    #load the LD matrix
    npz_file = '%s.npz'%(ld_prefix)
    try:
        R = load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s'%(npz_file))

    #create df_R and return it
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps

folder_ukbb = sys.argv[1]
chrom = int(sys.argv[2])
rsid_list = sys.argv[3].split(',')                                      
min_ld = int(sys.argv[4])

my_df_R, df_ld_snps = load_ld_npz(folder_ukbb)                                                              

print(f"matrix {os.path.dirname(folder_ukbb)} loaded")
print()

snps = df_ld_snps[df_ld_snps['SNP'].isin(rsid_list)]   
print(snps)
print()
snps = pd.DataFrame({'snp' : snps.index.tolist(), 'rsid' : snps['SNP'].tolist()})                                                               
print(snps)

high_ld_snps = []                  
for snp in snps['snp']:                 
    df_filtr = my_df_R.loc[snp, my_df_R.loc[snp] > (min_ld/100)]
    high_ld_snps.append(pd.DataFrame(df_filtr).T)
 
for df, rsid in zip(high_ld_snps, snps['rsid']):                            
    snp_list = df.columns.tolist()
    final_list = df_ld_snps[df_ld_snps.index.isin(snp_list)]
    print(f"snps list for {rsid} found")
    os.makedirs(f"../data/snp/highest-ld-{min_ld}/chrom-{chrom}", exist_ok=True)
    df.to_csv(f"../data/snp/highest-ld-{min_ld}/chrom-{chrom}/{rsid}-matrix.csv")
    final_list.to_csv(f"../data/snp/highest-ld-{min_ld}/chrom-{chrom}/{rsid}.csv", index=False)

print()
