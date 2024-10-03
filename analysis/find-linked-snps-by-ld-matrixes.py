import sys
import json
import csv
import pandas
import scipy.sparse as scipy_sparse

ld_matrixes_directory = sys.argv[1]
ld_value = float(sys.argv[2])

results_file_path = '../data/linked_snps.csv'

grouped_rsids_file_path = '../data/rsids-grouped-by-ld-prefixes.json'
rsids_grouped_by_ld_prefixes = []

with open(grouped_rsids_file_path, 'r') as grouped_rsids_file:
    file_content = grouped_rsids_file.read()
    rsids_grouped_by_ld_prefixes = json.loads(file_content)

def load_ld_matrix_by_ld_prefix(ld_prefix):
    snps_description_file_path = f'{ld_matrixes_directory}/{ld_prefix}.gz'
    snps_description = pandas.read_table(snps_description_file_path, sep='\s+')

    snps_description_rsids = snps_description['rsid'].astype(str)
    snps_description_rsids_and_bp = snps_description['rsid'].astype(str) + '.' + snps_description['position'].astype(str)

    ld_matrix_file_path = f'{ld_matrixes_directory}/{ld_prefix}.npz'

    try:
        ld_matrix_data = scipy_sparse.load_npz(ld_matrix_file_path).toarray()
        ld_matrix_data += ld_matrix_data.T
    except ValueError:
        raise IOError(f'Corrupt file: {ld_matrix_file_path}')
    
    ld_matrix = pandas.DataFrame(ld_matrix_data, index=snps_description_rsids_and_bp, columns=snps_description_rsids)

    return ld_matrix

def find_linked_snps_by_ld_prefix_and_rsid_list(ld_prefix, rsid_list):
    ld_matrix = load_ld_matrix_by_ld_prefix(ld_prefix)
    linked_snps = []

    for column_name, column_data in ld_matrix.items():
        if not column_name in rsid_list:
            continue

        for row_name, value in column_data.items():
            if float(value) <= ld_value:
                continue

            rs_id = row_name.split('.')[0]

            if rs_id == column_name:
                continue

            if not rs_id.startswith('rs'):
                continue

            chr = int(ld_prefix.replace('chr', '').split('_')[0])
            bp = int(row_name.split('.')[1])

            linked_snps.append({
                'RS_ID': rs_id,
                'CHR': chr,
                'BP': bp,
                'BASE_SNP': column_name
            })
    return linked_snps

linked_snps = []

for rsids_group in rsids_grouped_by_ld_prefixes:
    ld_prefix = rsids_group['ld_prefix']
    rsids = rsids_group['rs_ids']

    linked_snps += find_linked_snps_by_ld_prefix_and_rsid_list(ld_prefix, rsids)

dataframe = pandas.DataFrame(linked_snps)
dataframe.to_csv(results_file_path, index=False)
