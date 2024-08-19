import os
import re
import sys
import pandas as pd
import subprocess


#find ld file by bp in chromosome
def find_ld_file(bp, chrom, folder) :
    for filename in os.listdir(folder):
        if filename.endswith('.gz'):
            match = re.match(r'chr(\d+)_(\d+)_(\d+).gz', filename)
            if match:
                if chrom == int(match.group(1)):
                    begin = int(match.group(2))
                    end = int(match.group(3))
                    if begin <= bp <= end:
                        return folder + '/' + os.path.splitext(filename)[0]
    return None

def process_group(group):
    file_name = group['file_name'].iloc[0]
    rsids = group.iloc[:, 0].tolist()
    return file_name, rsids


#folder with LD-matrixes from aws
folder = sys.argv[1]
ld_value = int(sys.argv[2])
associations_path = '../data/snp/highest-ld-75'

for chrom in range(1, 23):
    rsid_file = f'../data/snp/by-chromosome/chromosome{chrom}.csv'
    existing_rsid_list = []
    df = pd.read_csv(rsid_file, names=['rsid', 'bp'])
    for index, row in df.iterrows():
        rsid = row['rsid']
        bp = row['bp']
        matrix_data_path = f'chrom-{chrom}/{rsid}-matrix.csv'

        if os.path.exists(os.path.join(associations_path, matrix_data_path)):
           # print(f"File for rsid={rsid} in chromosome {chrom} already exists.")
            existing_rsid_list.append(rsid)  
        else:
            print(f"File for rsid={rsid} in chromosome {chrom} does not exist.")

    df = df[~df['rsid'].isin(existing_rsid_list)]

    df['file_name'] = df.apply(lambda row: find_ld_file(row['bp'], chrom, folder), axis=1)

    grouped = df.groupby('file_name')

    for file_name, group_rsid in grouped:
        file_name, rsids = process_group(group_rsid)
        subprocess.run(f"python3 find-snps-with-ld.py {file_name} {chrom} {','.join(rsids)} {ld_value}", shell=True)
        
