import csv
import json
import os
import re
import requests
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
lock = threading.Lock()

def search_gwas_catalog(chromosome, main_rsid, list_rsids, not_found_rsids):
    file_path = f'../data/snp/associations/{main_rsid}.csv'
    if os.path.exists(file_path):
        print(f"File for {main_rsid} already exists. Skip")
        return
    all_data = []
    print(f'{main_rsid}: {list_rsids}')
    for rsid in list_rsids:
        print()
        url = f'https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chromosome}/associations/{rsid}'#?start=0&size=50'
        while url:
            print(url)
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                if '_embedded' not in data:
                    if 'p_value' in data and 'study_accession' in data and 'trait' in data:
                        trait = data['trait']
                        p_value = data['p_value']
                        study_accession = data['study_accession']
                        all_data.append({'chromosome': chromosome, 'rsid': rsid, 'trait': trait, 'study_accession': study_accession, 'p_value': p_value})
                    break
                associations = data['_embedded']['associations'].values()
                for association in associations:
                    trait = association['trait']
                    p_value = association['p_value']
                    study_accession = association['study_accession']
                    all_data.append({'chromosome': chromosome, 'rsid': rsid, 'trait': trait, 'study_accession': study_accession, 'p_value': p_value})
                links = data.get('_links')
                url = links.get('next', {}).get('href') if links else None
            else:
                url = None

    if all_data:
        if not os.path.exists('../data/snp/associations'):
            os.makedirs('../data/snp/associations')
        with open(file_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['chromosome', 'rsid', 'trait', 'study_accession', 'p_value'])
            writer.writeheader()
            for row in all_data:
                if row['p_value'] is not None and row['p_value'] <= 0.05 and row['p_value'] > 0:
                    writer.writerow(row)
        print(f'[INFO] Data for {main_rsid} in chromosome {chromosome} saved.')
    else:
        with lock:
            not_found_rsids.append(main_rsid) 
        print(f'[INFO] No data for {main_rsid} in chromosome {chromosome} in GWAS catalog.')
    return

def process_files_in_directory(directory):
    pattern = re.compile(r'rs\d+\.csv')
    futures = []
    not_found_rsids = []
    with ThreadPoolExecutor(max_workers=1) as executor:
        all_rsids = 0
        for i in range(1, 23):
            chromosome_folder = f'chrom-{i}'
            number_of_files = 0
            folder_path = os.path.join(directory, chromosome_folder)
            for root, dirs, files in os.walk(folder_path):
                for file in files:
                    if pattern.match(file):  
                        main_rsid = file.split('.')[0]
                        file_path = os.path.join(root, file)
                        list_rsids = []
                        number_of_files = number_of_files + 1
                        print(number_of_files)
                        with open(file_path, 'r') as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                rsid = row['SNP'] 
                                list_rsids.append(rsid)
                        search_gwas_catalog(i, main_rsid, list_rsids, not_found_rsids)
                        futures.append(executor.submit(search_gwas_catalog, i, main_rsid, list_rsids, not_found_rsids))
            print(f'[INFO] for chrom {i} there is {number_of_files} rsids')
            all_rsids = all_rsids + number_of_files
        print(f'[INFO] {all_rsids} rsids')
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"An error occurred: {e}")
    not_found_file_path = '../data/snp/associations/not-found-rdids.csv'
    with open(not_found_file_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['rsid'])
            for rsid in not_found_rsids:
                writer.writerow([rsid])
    print(f"[INFO] RSIDs not found saved to {not_found_file_path}")

directory_path = '../data/snp/highest-ld-75'
process_files_in_directory(directory_path)
print(f'\n all rsids done\n')
