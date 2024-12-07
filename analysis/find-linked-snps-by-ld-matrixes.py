import sys
import json
import pandas
import scipy.sparse as scipy_sparse

ld_matrixes_directory = sys.argv[1]
ld_value = float(sys.argv[2])

results_file_path = '../data/snps-found-via-ld-matrixes.csv'
grouped_rsids_file_path = '../data/rsids-grouped-by-ld-prefixes.json'
completed_rsids_groups_file_path = '../data/completed-rsids-groups.json'

def read_json_file(file_path):
    try:
        with open(file_path, 'r') as file:
            file_content = file.read()
            try:
                return json.loads(file_content)
            except json.decoder.JSONDecodeError:
                return []
    except FileNotFoundError:
        return []

def array_of_dictionaries_from_csv_file(file_path):
    try:
        snps_file = pandas.read_csv(file_path)
    except (FileNotFoundError, pandas.errors.EmptyDataError):
        return []
    
    return [dict(snp[1]) for snp in snps_file.iterrows()]

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

def find_linked_snps_by_ld_prefix_and_rsid_list(ld_prefix, rsid_list, already_added_snps):
    ld_matrix = load_ld_matrix_by_ld_prefix(ld_prefix)
    linked_snps = []

    for column_name, column_data in ld_matrix.items():
        if not column_name in rsid_list:
            continue

        for row_name, value in column_data.items():
            if float(value) <= ld_value:
                continue

            rs_id = row_name.split('.')[0]

            if not rs_id.startswith('rs'):
                continue

            chr = int(ld_prefix.replace('chr', '').split('_')[0])
            bp = int(row_name.split('.')[1])

            if any([
                rs_id == already_added_snp['RS_ID'] and column_name == already_added_snp['BASE_SNP']
                for already_added_snp in already_added_snps
                ]):
                continue

            linked_snps.append({
                'RS_ID': rs_id,
                'CHR': chr,
                'BP': bp,
                'BASE_SNP': column_name,
                'LD_VALUE': float(value)
            })
    return linked_snps

found_snps = array_of_dictionaries_from_csv_file(results_file_path)
rsids_grouped_by_ld_prefixes = read_json_file(grouped_rsids_file_path)
completed_rsids_groups = read_json_file(completed_rsids_groups_file_path)

for rsids_group in rsids_grouped_by_ld_prefixes:
    ld_prefix = rsids_group['ld_prefix']
    rsids = rsids_group['rs_ids']

    if (ld_prefix in completed_rsids_groups):
        continue
 
    found_snps += find_linked_snps_by_ld_prefix_and_rsid_list(ld_prefix, rsids, found_snps)
    
    dataframe = pandas.DataFrame(found_snps)
    dataframe.to_csv(results_file_path, index=False)

    completed_rsids_groups.append(ld_prefix)
    with open(completed_rsids_groups_file_path, 'w') as file:
        json.dump(completed_rsids_groups, file, indent=4)
