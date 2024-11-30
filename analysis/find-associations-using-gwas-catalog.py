import sys
import csv
import requests
import re
import pandas

p_value = float(sys.argv[1])
snps_file_path = '../data/snps-found-via-ld-matrixes.csv'
results_file_path = '../data/associations-found-by-gwas-catalog.csv'
snps_list = []

def get_api_response(chr, rs_id, start_query_parameter = 0, size_query_parameter = 20):
    query_arguments = f'start={start_query_parameter}&size={size_query_parameter}&p_lower=0&p_upper={p_value}'
    api_url = f'https://www.ebi.ac.uk/gwas/summary-statistics/api/chromosomes/{chr}/associations/{rs_id}?{query_arguments}'

    api_response = requests.get(api_url, verify=True)

    if (api_response.status_code == 404):
        return 'Not found'

    if 'next' in api_response.json()['_links'] and '_embedded' in api_response.json()['_links']:
        api_response_next_url = api_response.json()['_links']['next']['href']
        api_response_next_url_start_value = re.search(r'start=([0-9]|\.)+', api_response_next_url).group(0).split('=')[1]
    else:
        api_response_next_url_start_value = None

    if '_embedded' in api_response.json():
        embedded_associations = api_response.json()['_embedded']['associations']
    else:
        embedded_associations = {
            '0': api_response.json()
        }

    api_response_formated = {
        "associations": [
            embedded_associations
        ],
        "next_page_data": {
            "start": api_response_next_url_start_value
        }
    }

    return api_response_formated

def get_associations_by_chr_and_rsid(chr, rsid):
    start = 0
    size = 20

    whole_response = []
    while True:
        api_response = get_api_response(chr, rsid, start, size)

        if (api_response == 'Not found'):
            return

        for association in api_response['associations'][0].values():
            whole_response.append(association)

        if (api_response['next_page_data']['start'] == None):
            break

        start = api_response['next_page_data']['start']

    return whole_response

def read_snps_csv_and_write_data_to_snps_list():
    snps_csv = open(snps_file_path, 'r', newline='')

    csv_reader = csv.reader(snps_csv)

    csv_header = next(csv_reader)
    csv_header_formated = [value.strip() for value in csv_header]

    rs_id_index = csv_header_formated.index('RS_ID')
    chr_index = csv_header_formated.index('CHR')
    bp_index = csv_header_formated.index('BP')
    base_snp_index = csv_header_formated.index('BASE_SNP')
    ld_value_index = csv_header_formated.index('LD_VALUE')

    for row in csv_reader:
        row_formated = [value.strip() for value in row]

        snp_object = {
            'rs_id': row_formated[rs_id_index],
            'chr':  int(row_formated[chr_index]),
            'bp':   int(row_formated[bp_index]),
            'base_snp': row_formated[base_snp_index],
            'ld_value': row_formated[ld_value_index]
        }

        snps_list.append(snp_object)
               
    snps_csv.close()

read_snps_csv_and_write_data_to_snps_list()

required_columns = ['RS_ID', 'CHR', 'BP', 'BASE_SNP', 'LD_VALUE', 'P_VALUE', 'BETA', 'ODDS_RATIO', 'EFO_TRAIT', 'STUDY_ID', 'API_NAME']

try:
    results_file_dataframe = pandas.read_csv(results_file_path)
except (pandas.errors.EmptyDataError, FileNotFoundError) as error:
    with open(results_file_path, 'w', newline='') as f:
        f.write(','.join(required_columns))
    results_file_dataframe = pandas.read_csv(results_file_path)

for required_column in required_columns:
    if required_column in results_file_dataframe.columns:
        continue
    results_file_dataframe.insert(len(results_file_dataframe.columns), required_column, None)
    
results_file_dataframe.to_csv(results_file_path, index=False)

for snp in snps_list:
    rs_id = snp['rs_id']
    chr = snp['chr']
    bp = snp['bp']
    base_snp = snp['base_snp']
    ld_value = snp['ld_value']

    associations = get_associations_by_chr_and_rsid(chr, rs_id)
    if associations == None:
        continue

    association_rows = []
    for association in associations:
        association_rows.append({
            'RS_ID': rs_id,
            'CHR': chr,
            'BP': bp,
            'BASE_SNP': base_snp,
            'LD_VALUE': ld_value,
            'P_VALUE': float(association['p_value']),
            'BETA': float(association['beta']) if association['beta'] is not None else 'null',
            'ODDS_RATIO': float(association['odds_ratio']) if association['odds_ratio'] is not None else 'null',
            'EFO_TRAIT': association['trait'][0],
            'STUDY_ID': association['study_accession'],
            'API_NAME': 'gwas-catalog'
        })
    results_file_dataframe = pandas.concat([results_file_dataframe, pandas.DataFrame(association_rows)]).reset_index(drop=True)
    results_file_dataframe.to_csv(results_file_path, index=False)

# TODO: connect open targets
# TODO: connect IEU Open GWAS project 
# TODO: create script that lookups the EFO of the snps and prints the traits short description