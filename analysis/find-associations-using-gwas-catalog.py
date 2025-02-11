import sys
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

    print(api_url)
    api_response = requests.get(api_url, verify=True)
    print(api_response)

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

def go_through_snps_and_save_associations_to_file(snps_file_path, results_file_path):
    snps_file = pandas.read_csv(snps_file_path)
    try:
        results_file = pandas.read_csv(results_file_path)
    except (FileNotFoundError, pandas.errors.EmptyDataError):
        results_file = pandas.DataFrame(columns=['RS_ID', 'CHR', 'BP', 'BASE_SNP', 'LD_VALUE', 'P_VALUE', 'BETA', 'ODDS_RATIO', 'EFO_TRAIT', 'STUDY_ID'])
    
    for snp in snps_file.iterrows():
        snp_dict = dict(snp[1])

        if snp_dict['RS_ID'] in results_file['RS_ID'].unique():
            continue

        snp_associations = get_associations_by_chr_and_rsid(snp_dict['CHR'], snp_dict['RS_ID'])
        if snp_associations == None:
            continue

        snp_associations_formated = []
        for snp_association in snp_associations:
            for trait in snp_association['trait']:
                snp_associations_formated.append({
                    'RS_ID': snp_dict['RS_ID'],
                    'CHR': snp_dict['CHR'],
                    'BP': snp_dict['BP'],
                    'BASE_SNP': snp_dict['BASE_SNP'],
                    'LD_VALUE': snp_dict['LD_VALUE'],
                    'P_VALUE': float(snp_association['p_value']),
                    'BETA': float(snp_association['beta']) if snp_association['beta'] is not None else 'null',
                    'ODDS_RATIO': float(snp_association['odds_ratio']) if snp_association['odds_ratio'] is not None else 'null',
                    'EFO_TRAIT': trait,
                    'STUDY_ID': snp_association['study_accession']
                })
    
        results_file = pandas.concat([results_file, pandas.DataFrame(snp_associations_formated)]).reset_index(drop=True)
        results_file.to_csv(results_file_path, index=False)

go_through_snps_and_save_associations_to_file(snps_file_path, results_file_path)