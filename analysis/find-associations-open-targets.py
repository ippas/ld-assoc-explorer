import csv
import json
import os
import re
import requests
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
lock = threading.Lock()

base_url = "https://api.genetics.opentargets.org/graphql"

search_query = """
  query searchRsId($queryString: String!) {
    search(queryString: $queryString) {
      variants {
        id
      }
    }
  }
"""

variant_info_query = """
  query variantInfo($variantId: String!) {
    indexVariantsAndStudiesForTagVariant(variantId: $variantId) {
      associations {
        study {
          studyId
          traitReported
          traitCategory
          pmid
          pubAuthor
          pubDate
        }
        indexVariant {
          id
          rsId
          refAllele
          altAllele
          nearestGeneDistance
          nearestGene {
            symbol
          }
        }
        pval
        nTotal
        nCases
        overallR2
        afr1000GProp
        amr1000GProp
        eas1000GProp
        eur1000GProp
        sas1000GProp
        log10Abf
        posteriorProbability
        pvalMantissa
        pvalExponent
        oddsRatio
        oddsRatioCILower
        oddsRatioCIUpper
        beta
        direction
        betaCILower
        betaCIUpper
      }
    }
  }
"""


def get_variant_id(rsid, retries=5, backoff_factor=1):
    variables = {"queryString": rsid}
    for attempt in range(retries):
        try:
            response = requests.post(base_url, json={"query": search_query, "variables": variables})
            if response.status_code == 200:
                response_json = response.json()
                variants = response_json.get("data", {}).get("search", {}).get("variants", [])
                if variants:
                    return variants[0].get("id")
                return None
            else:
                raise Exception(f"Query failed with status code {response.status_code}, response: {response.text}")
        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(backoff_factor * (2 ** attempt))
    raise Exception(f"All {retries} attempts failed for rsid: {rsid}")

def get_variant_associations(variant_id, retries=5, backoff_factor=1):
    variables = {"variantId": variant_id}
    for attempt in range(retries):
        try:
            response = requests.post(base_url, json={"query": variant_info_query, "variables": variables})
            if response.status_code == 200:
                return response.json()
            else:
                raise Exception(f"Query failed with status code {response.status_code}, response: {response.text}")
        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(backoff_factor * (2 ** attempt))
    raise Exception(f"All {retries} attempts failed for variant ID: {variant_id}")


def search_open_targets(chromosome, main_rsid, list_rsids, not_found_rsids):
    file_path = f'../data/snp/associations-open-targets/{main_rsid}.csv'
    
    if os.path.exists(file_path):
        print(f"File for {main_rsid} already exists. Skip")
        return
    
    all_data = []
    print(f'Processing {main_rsid}: {list_rsids}')
    
    for rsid in list_rsids:
        try:
            if rsid.startswith("rs"):
                variant_id = get_variant_id(rsid)
            else:
                variant_id = rsid
                variant_id = variant_id.replace(':', '_')
            print(variant_id) 
            if variant_id:
                info = get_variant_associations(variant_id)
                associations = info.get("data", {}).get("indexVariantsAndStudiesForTagVariant", {}).get("associations", [])
                
                for assoc in associations:
                    result = {
                        'chromosome': chromosome,
                        'rsid': rsid,
                        'trait': assoc["study"]["traitReported"],
                        'study_accession': assoc["study"]["studyId"],
                        'p_value': assoc["pval"]
                    }
                    all_data.append(result)
        except Exception as e:
            print(f"Error retrieving data for {rsid}: {e}")
            continue
    
    if all_data:
        if not os.path.exists('../data/snp/associations-open-targets'):
            os.makedirs('../data/snp/associations-open-targets')
        
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
            print(f'[INFO] No data for {main_rsid} in chromosome {chromosome} in OpenTargets.')

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
            
                        search_open_targets(i, main_rsid, list_rsids, not_found_rsids)
        
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





