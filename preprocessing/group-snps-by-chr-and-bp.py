import sys
import os
import pandas
import json

SNPS_FILE_PATH = sys.argv[1]

LD_MATRIXES_DIRECTORY = sys.argv[2]
LD_FILE_NAMES = os.listdir(LD_MATRIXES_DIRECTORY)

RESULTS_FILE_PATH = '../data/rsids-grouped-by-ld-prefixes.json'

def error(message: str) -> None:
    script_name = sys.argv[0].split('/')[-1]
    print(f"{script_name}: {message}", file=sys.stderr)
    sys.exit(1)

def ld_prefixes_from_ld_files(ld_file_names: list) -> list:
    ld_prefixes = []
    for ld_file in ld_file_names:
        if not ld_file.startswith('chr'):
            continue
        if not ld_file.endswith('.gz'):
            continue
        ld_prefixes.append(ld_file.replace('.gz', ''))
    return ld_prefixes

def snps_list_from_snps_file(snps_file_path: str) -> list:
    snps_file_name = snps_file_path.split('/')[-1]

    try:
        snps_file = pandas.read_csv(snps_file_path)
    except pandas.errors.EmptyDataError:
        error(f'{snps_file_name}: cannot be empty')

    snps_file_columns = snps_file.columns
    required_columns = ['RS_ID', 'CHR', 'BP']

    if any([required_column not in snps_file_columns for required_column in required_columns]):
        error(f'{snps_file_name}: must contain header with following columns: {required_columns}')
    
    return [dict(snp[1]) for snp in snps_file.iterrows()]

def matching_ld_prefixes_by_chr_and_bp(ld_prefixes: list, chr: int, bp: int) -> list:
    matching_ld_prefixes = []

    for ld_prefix in ld_prefixes:
        ld_prefix_formated = [int(value) for value in ld_prefix.replace('chr', '').split('_')]

        if ld_prefix_formated[0] != int(chr):
            continue
        if ld_prefix_formated[1] > int(bp):
            continue
        if ld_prefix_formated[2] < int(bp):
            continue

        matching_ld_prefixes.append(ld_prefix)

    return matching_ld_prefixes

def append_rsid_grouped_by_ld_prefix_to_list(grouped_rsids: list, ld_prefix: str, snp: str) -> list:
    new_group = {
        'ld_prefix': ld_prefix,
        'rs_ids': [
            snp
        ]
    }

    if len(grouped_rsids) == 0:
        grouped_rsids.append(new_group)
        return grouped_rsids
    
    for group in grouped_rsids:
        if group['ld_prefix'] != ld_prefix:
            continue

        if snp in group['rs_ids']:
            return grouped_rsids

        group['rs_ids'].append(snp)
        return grouped_rsids
        
    grouped_rsids.append(new_group)
    return grouped_rsids

def main() -> None:
    ld_prefixes = ld_prefixes_from_ld_files(LD_FILE_NAMES)
    snps_list = snps_list_from_snps_file(SNPS_FILE_PATH)

    grouped_rsids = []

    for snp in snps_list:
        matching_ld_prefixes = matching_ld_prefixes_by_chr_and_bp(ld_prefixes, snp['CHR'], snp['BP'])
        
        rsid = snp['RS_ID']
        for matching_ld_prefix in matching_ld_prefixes:
            grouped_rsids = append_rsid_grouped_by_ld_prefix_to_list(grouped_rsids, matching_ld_prefix, rsid)

    with open(RESULTS_FILE_PATH, 'w') as file:
        file.write(json.dumps(grouped_rsids, indent=4))

if __name__ == "__main__":
    main()