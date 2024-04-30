import json
import os
import pandas as pd
import csv

snp_table = pd.read_csv('../data/snp/nearest_with_rare_march.csv')
snp_table = snp_table.sort_values(by=['chr', 'bp']).reset_index(drop=True)

grouped = snp_table.groupby('chr')

for chr_num, group_data in grouped:
    chromosome_folder ='../data/snp/by-chromosome/' 
    file_path = os.path.join(chromosome_folder, f'chromosome{chr_num}.csv')
    with open(file_path, 'a', newline='') as f:
        file = csv.writer(f, delimiter=',')
        for _, row in group_data.iterrows():
            file.writerow([row['rsID'], row['bp']])
