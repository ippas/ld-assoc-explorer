import csv
import json
import os
import re
import requests


def search_gwas_catalog(chromosome, main_rsid, list_rsids):
    all_data = []
    print(f'{main_rsid}: {list_rsids}')
    for rsid in list_rsids:
        url = f'https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chromosome}/associations/{rsid}'
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            for association in data['_embedded']['associations'].values():
                trait = association['trait']
                p_value = association['p_value']
                study_accession = association['study_accession']
                all_data.append({'chromosome': chromosome, 'rsid': rsid, 'trait': trait, 'study_accession': study_accession, 'p_value': p_value})
        else:
            print(f'{rsid} not found')  
    print()
    if all_data:
        if not os.path.exists('../data/snp/associations'):
            os.makedirs('../data/snp/associations')
        file_path = f'../data/snp/associations/{main_rsid}.csv'
        with open(file_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['chromosome','rsid', 'trait', 'study_accession', 'p_value'])
            writer.writeheader()
            for row in all_data:
                if row['p_value'] != None:
                    if row['p_value'] <= 0.05:
                        writer.writerow(row)

def process_files_in_directory(directory):
    pattern = re.compile(r'rs\d+\.csv')
    for i in range(8, 9):
        chromosome_folder = f'chrom-{i}'
        folder_path = os.path.join(directory, chromosome_folder)
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                if pattern.match(file):  
                    main_rsid = file.split('.')[0]
                    file_path = os.path.join(root, file)
                    list_rsids = []
                    with open(file_path, 'r') as f:
                        reader = csv.DictReader(f)
                        for row in reader:
                            rsid = row['SNP'] 
                            list_rsids.append(rsid)
                    search_gwas_catalog(i, main_rsid, list_rsids)

directory_path = '../data/snp/highest-ld-75'
process_files_in_directory(directory_path)