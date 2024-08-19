import ieugwaspy as igd
import pandas as pd
import json
import csv
import os
import re
import requests
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

lock = threading.Lock()
ids = list(igd.gwasinfo().keys())

def search_opengwas(chrom, main_rsid, rsids, not_found_rsids):
    list_size = len(ids)
    chunk_size = 500    
    file_path = f'/app/data/snp/associations-opengwas/{main_rsid}.csv'
    results = []

    if os.path.exists(file_path):
        print(f"File for {main_rsid} already exists. Skip")
        return

    for i in range(0, list_size, chunk_size):
        chunk = ids[i:i + chunk_size]
        assoc = igd.associations(id = chunk, variant = rsids)
        if assoc:
            print(f'[INFO] associations was found')
            for value in assoc:
                pvalue = value.get("p")
            
                if pvalue is not None and pvalue < 0.05 and pvalue > 0:
                    result = {
                        "chromosome": chrom,
                        "rsid": value.get("rsid"),
                        "trait": value.get("trait"),
                        "study_accession": value.get("id"),
                        "p_value": pvalue
                    }
                    results.append(result)
        
    if results:
        if not os.path.exists('/app/data/snp/associations-opengwas'):
            os.makedirs('/app/data/snp/associations-opengwas')
        with open(file_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['chromosome', 'rsid', 'trait', 'study_accession', 'p_value'])
            writer.writeheader()
            for row in results:
                writer.writerow(row)
        print(f'[INFO] Data for {main_rsid} in chromosome {chromosome} saved.')
    else:
        with lock:
            not_found_rsids.append(main_rsid)
            print(f'[info] Data for {main_rsid} in {chromosome} not found.')
    
def process_files_in_directory(directory):
    pattern = re.compile(r'rs\d+\.csv')
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
            
                        search_opengwas(i, main_rsid, list_rsids, not_found_rsids)
        
            print(f'[INFO] for chrom {i} there is {number_of_files} rsids')
            all_rsids = all_rsids + number_of_files
        
        print(f'[INFO] {all_rsids} rsids')
        
    
    not_found_file_path = '/app/data/snp/associations/not-found-rdids.csv'
    with open(not_found_file_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['rsid'])
            for rsid in not_found_rsids:
                writer.writerow([rsid])
    print(f"[INFO] RSIDs not found saved to {not_found_file_path}")

directory_path = '/app/data/snp/highest-ld-75'
process_files_in_directory(directory_path)
print(f'\n all rsids done\n')
