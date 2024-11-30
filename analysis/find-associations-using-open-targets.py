import sys
import requests
import csv
import pandas


ARG_P_VALUE = float(sys.argv[1])
snps_file_path = '../data/snps-found-via-ld-matrixes.csv'
results_file_path = '../data/associations-found-by-open-targets.csv'

API_URL = "https://api.genetics.opentargets.org/graphql"

def get_variant_id_from_rs_id(rs_id: str) -> str:
    search_query = """
	query searchRsId($rsId: String!) {
		search(queryString: $rsId) {
	    variants {
	      id
	    }
	  }
	}
    """
    variables = {"rsId": rs_id}

    search_result = requests.post(API_URL, json={'query': search_query, 'variables':variables}).json()

    try:
        variant_id = search_result['data']['search']['variants'][0]['id']
    except IndexError:
        variant_id = None

    return variant_id

def get_variant_data(variant_id: str) -> dict:
    variant_query = """
    query PheWASQuery($variantId: String!) {
    pheWAS(variantId: $variantId) {
        totalGWASStudies
        associations {
        study {
            studyId
            pmid
            source
            traitEfos
            __typename
        }
        pval
        beta
        oddsRatio
        __typename
        }
        __typename
    }
    }
    """
    variables = {"variantId": variant_id}

    variant_result = requests.post(API_URL, json={'query': variant_query, 'variables': variables}).json()
    return variant_result

def get_rs_id_associations_filtered_by_p_value(rs_id: str, p_lower: float, p_upper: float) -> list:
    variant_id = get_variant_id_from_rs_id(rs_id)

    if (variant_id == None):
        return []

    variant_data = get_variant_data(variant_id)

    rs_id_associations = []

    variant_associations = variant_data['data']['pheWAS']['associations']
    for association in variant_associations:
        if (association['pval'] < p_lower or association['pval'] > p_upper):
            continue

        rs_id_associations.append({
            'study_id': association['study']['studyId'],
            'pmid': association['study']['pmid'],
            'efos': association['study']['traitEfos'],
            'p_value': float(association['pval']) if (association['pval']) is not None else 'null',
            'beta': float(association['beta']) if (association['beta']) is not None else 'null',
            'odds_ratio': float(association['oddsRatio']) if (association['oddsRatio']) is not None else 'null'
        })

    return rs_id_associations

def get_rs_ids_from_snps_csv() -> list:
    snps_csv = open(snps_file_path, 'r', newline='')

    csv_reader = csv.reader(snps_csv)

    csv_header = next(csv_reader)
    csv_header_formated = [value.strip() for value in csv_header]

    rs_id_index = csv_header_formated.index('RS_ID')
    chr_index = csv_header_formated.index('CHR')
    bp_index = csv_header_formated.index('BP')
    base_snp_index = csv_header_formated.index('BASE_SNP')
    ld_value_index = csv_header_formated.index('LD_VALUE')

    snps_list = []

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
    
    return snps_list

required_columns = ['RS_ID', 'CHR', 'BP', 'BASE_SNP', 'LD_VALUE', 'P_VALUE', 'BETA', 'ODDS_RATIO', 'EFO_TRAIT', 'STUDY_ID', 'PMID', 'API_NAME']

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

snps = get_rs_ids_from_snps_csv()
for snp in snps:
    rs_id_associations = get_rs_id_associations_filtered_by_p_value(snp['rs_id'], 0, ARG_P_VALUE)
    formated_associations = []
    for association in rs_id_associations:
        for index, efo in enumerate(association['efos']):
            formated_associations.append({
                'RS_ID': snp['rs_id'],
                'CHR': snp['chr'],
                'BP': snp['bp'],
                'BASE_SNP': snp['base_snp'],
                'LD_VALUE': snp['ld_value'],
                'P_VALUE': float(association['p_value']),
                'BETA': float(association['beta']) if association['beta'] != 'null' else 'null',
                'ODDS_RATIO': float(association['odds_ratio']) if association['odds_ratio'] != 'null' else 'null',
                'EFO_TRAIT': efo,
                'STUDY_ID': association['study_id'],
                'PMID': association['pmid'],
                'API_NAME': 'open-targets'
            })
    results_file_dataframe = pandas.concat([results_file_dataframe, pandas.DataFrame(formated_associations)]).reset_index(drop=True)
    results_file_dataframe.to_csv(results_file_path, index=False)