import sys
import requests
import csv
import pandas

ARG_P_VALUE = float(sys.argv[1])
FOUND_SNPS_FILE_PATH = '../data/snps-found-via-ld-matrixes.csv'
RESULTS_FILE_PATH = '../data/associations-found-by-open-targets.csv'
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

    for _ in range(5):
        search_result = requests.post(API_URL, json={'query': search_query, 'variables':variables})

        if search_result.status_code != 200:
            continue

        try:
            return search_result.json()['data']['search']['variants'][0]['id']
        except:
            return None

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

    for _ in range(5):
        variant_result = requests.post(API_URL, json={'query': variant_query, 'variables': variables})

        if variant_result.status_code != 200:
            continue

        try:
            return variant_result.json()
        except:
            return None

def get_rs_id_associations_filtered_by_p_value(rs_id: str, p_lower: float, p_upper: float) -> list:
    variant_id = get_variant_id_from_rs_id(rs_id)
    if (variant_id == None):
        return []

    variant_data = get_variant_data(variant_id)
    if (variant_data == None):
        return []

    rs_id_associations = []

    variant_associations = variant_data['data']['pheWAS']['associations']
    for association in variant_associations:
        if (float(association['pval']) < p_lower or float(association['pval']) > p_upper):
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

def array_of_dictionaries_from_csv_file(file_path: str) -> list:
    try:
        snps_file = pandas.read_csv(file_path)
    except (FileNotFoundError, pandas.errors.EmptyDataError):
        return []
    
    return [dict(snp[1]) for snp in snps_file.iterrows()]

def save_progress(found_associations: list) -> None:
    dataframe = pandas.DataFrame(found_associations)
    dataframe.to_csv(RESULTS_FILE_PATH, index=False)

def main() -> None:
    snps_found_via_ld_matrixes = array_of_dictionaries_from_csv_file(FOUND_SNPS_FILE_PATH)
    found_associations = array_of_dictionaries_from_csv_file(RESULTS_FILE_PATH)

    completed_snps = [
        (snp['RS_ID'], snp['BASE_SNP'])
        for snp in found_associations
    ]

    snps_to_process = [
        snp
        for snp in snps_found_via_ld_matrixes
        if (snp['RS_ID'], snp['BASE_SNP']) not in completed_snps
    ]

    for snp in snps_to_process:
        rs_id_associations = get_rs_id_associations_filtered_by_p_value(snp['RS_ID'], 0, ARG_P_VALUE)
        if rs_id_associations == []:
            found_associations.append({
                    'RS_ID': snp['RS_ID'],
                    'CHR': snp['CHR'],
                    'BP': snp['BP'],
                    'BASE_SNP': snp['BASE_SNP'],
                    'LD_VALUE': snp['LD_VALUE']
                })
            save_progress(found_associations)
            continue

        for association in rs_id_associations:
            for efo_index in range(max(1,len(association['efos']))):
                found_associations.append({
                    'RS_ID': snp['RS_ID'],
                    'CHR': snp['CHR'],
                    'BP': snp['BP'],
                    'BASE_SNP': snp['BASE_SNP'],
                    'LD_VALUE': snp['LD_VALUE'],
                    'P_VALUE': float(association['p_value']),
                    'BETA': float(association['beta']) if association['beta'] != 'null' else 'null',
                    'ODDS_RATIO': float(association['odds_ratio']) if association['odds_ratio'] != 'null' else 'null',
                    'EFO_TRAIT': association['efos'][efo_index] if len(association['efos']) > 0 else 'null',
                    'STUDY_ID': association['study_id'],
                    'PMID': association['pmid']
                })

        save_progress(found_associations)

if __name__ == "__main__":
    main()