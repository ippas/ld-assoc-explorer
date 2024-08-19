import os
import sys
import pandas as pd
import csv

filepath = sys.argv[1]
snp_table = pd.read_csv(filepath)

rsids = snp_table['rsID']
chrs = snp_table['chr']
bps = snp_table['bp']
i = 0

os.makedirs('../data/snp', exist_ok=True)
os.makedirs(f'../data/snp/by-chromosome', exist_ok=True)

for i in range(len(rsids)):
    rsid = rsids[i]
    chr_num = chrs[i]
    bp = bps[i]
    chromosome_folder ='../data/snp/by-chromosome/' 
    with open(os.path.join(chromosome_folder, f'chromosome{chr_num}.csv'), 'a', newline='') as f:
        file = csv.writer(f, delimiter=',')
        file.writerow([rsid, bp])
