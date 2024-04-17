import json
import pandas as pd
import requests

snp_table = pd.read_csv('data/snp/nearest_with_rare_march.csv')

rsids = snp_table['rsID']

def search_gwas_catalog(rsid):
    url = f'https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/associations'
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        with open(f'data/snp/associations/{rsid}_gwas_data.json', 'w') as f:
            json.dump(data, f, indent=4)

for rsid in rsids:
    search_gwas_catalog(rsid)
