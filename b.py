import json
import os
import pandas as pd
import csv

snp_table = pd.read_csv('data/snp/nearest_with_rare_march.csv')

rsids = snp_table['rsID']
chrs = snp_table['chr']
i = 0
for rsid in rsids:
    with open(f'data/snp/by_chromosome/chromosome{chrs[i]}.csv', 'w') as f:
        file = csv.writer(f, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        file.writerow([f'{rsid}'])
    i = i + 1
        