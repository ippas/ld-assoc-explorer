import csv
from itertools import groupby
import os
import re
import pandas as pd
import requests

def find_traits(reader):
    answer = []
    for row in reader:
        trait = row['trait'][2:-2]
        url= f'https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{trait}'
        response = requests.get(url)
        if response.status_code == 200:
            response = response.json()
            answer.append(response['trait'])
        else: 
            answer.append(trait)
    return answer

def process_files_in_directory(source_directory, answer_file):
    pattern = re.compile(r'rs\d+\.csv')
    with open(answer_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['rsid', 'traits'])
        writer.writeheader()
        for root, dirs, files in os.walk(source_directory):
            for file in files:
                print(file)
                if pattern.match(file):  
                    main_rsid = file.split('.')[0]
                    file_path = os.path.join(source_directory, file)
                    traits_names = []
                    with open(file_path, 'r') as f:
                        reader = csv.DictReader(f)
                        traits_names = find_traits(reader)
                    if len(traits_names) != 0:
                        traits_names = list(set(traits_names))
                        print(traits_names)
                        writer.writerow({'rsid': main_rsid, 'traits': ', '.join(traits_names)})
                print()
    f.close()

directory = '../data/snp/associations'
result = '../results/traits.csv'
no_results = '../results/no-traits.csv'

process_files_in_directory(directory, result)

rsid1 = pd.read_csv('../data/snp/nearest_with_rare_march.csv')['rsID']
rsid2 = pd.read_csv(result)['rsid']
diff = list(set(rsid1) - set(rsid2))

fieldnames = ['rsid']
with open(no_results, mode='w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in diff:
        writer.writerow({'rsid': row})
