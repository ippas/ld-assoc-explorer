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
for chrom in range(1, 23):
    rsid_file = f'../data/snp/by-chromosome/chromosome{chrom}.csv'

    df = pd.read_csv(rsid_file, names=['rsid', 'bp'])

    df['file_name'] = df.apply(lambda row: find_ld_file(row['bp'], chrom, folder), axis=1)

    grouped = df.groupby('file_name')

    for file_name, group_rsid in grouped:
        file_name, rsids = process_group(group_rsid)
        subprocess.run(f"python3 find-snps-with-ld.py {file_name} {chrom} {','.join(rsids)}", shell=True)
        