import sys
import os
import csv
import json

snps_file_path = sys.argv[1]
ld_matrixes_directory = sys.argv[2]

results_file_path = '../data/rsids-grouped-by-ld-prefixes.json'

snps_file_name = snps_file_path.split('/')[-1]
ld_files = os.listdir(ld_matrixes_directory)

ld_prefixes = []
snps_list = []
rsids_grouped_by_ld_prefixes= []

def error(message):
    script_name = sys.argv[0].split('/')[-1]
    print(f"{script_name}: {message}", file=sys.stderr)
    sys.exit(1)

def ld_files_to_ld_prefixes():
    for ld_file in ld_files:
        if not ld_file.startswith('chr'):
            continue
        if not ld_file.endswith('.gz'):
            continue
        ld_prefixes.append(ld_file.replace('.gz', ''))

def read_snps_csv_and_write_data_to_snps_list():
    snps_csv = open(snps_file_path, 'r', newline='')

    csv_reader = csv.reader(snps_csv)

    csv_header = next(csv_reader)
    csv_header_formated = [value.strip() for value in csv_header]
    needed_columns = ['RS_ID', 'CHR', 'BP']

    if not all(column in csv_header_formated for column in needed_columns):
        error(f'{snps_file_name}: must contain header with following columns: {needed_columns}')

    rs_id_index = csv_header_formated.index('RS_ID')
    chr_index = csv_header_formated.index('CHR')
    bp_index = csv_header_formated.index('BP')

    for row in csv_reader:
        row_formated = [value.strip() for value in row]

        rs_id = row_formated[rs_id_index]

        if not rs_id.startswith('rs'):
            rs_id = f'rs{rs_id}'

        try:
            snp_object = {
                'rs_id': rs_id,
                'chr':  int(row_formated[chr_index]),
                'bp':   int(row_formated[bp_index])
            }
        except ValueError as e:
            error(f'{snps_file_name}: {e}')

        snps_list.append(snp_object)
               
    snps_csv.close()

def ld_prefixes_by_chr_and_bp(chr, bp):
    founded_ld_prefixes = []

    for ld_prefix in ld_prefixes:
        ld_prefix_formated = [int(value) for value in ld_prefix.replace('chr', '').split('_')]

        if ld_prefix_formated[0] != chr:
            continue
        if ld_prefix_formated[1] > bp:
            continue
        if ld_prefix_formated[2] < bp:
            continue

        founded_ld_prefixes.append(ld_prefix)

    return founded_ld_prefixes

def add_snp_grouped_by_ld_prefix(ld_prefix, snp):
    new_object = {
        'ld_prefix': ld_prefix,
        'rs_ids': [
            snp
        ]
    }

    if len(rsids_grouped_by_ld_prefixes) == 0:
        rsids_grouped_by_ld_prefixes.append(new_object)
        return

    for object in rsids_grouped_by_ld_prefixes:
        if object['ld_prefix'] != ld_prefix:
            continue

        object['rs_ids'].append(snp)
        return

    rsids_grouped_by_ld_prefixes.append(new_object)

ld_files_to_ld_prefixes()
read_snps_csv_and_write_data_to_snps_list()

for snp in snps_list:
    for ld_prefix in ld_prefixes_by_chr_and_bp(snp['chr'], snp['bp']):
        add_snp_grouped_by_ld_prefix(ld_prefix, snp['rs_id'])

with open(results_file_path, 'w') as file:
    file.write(json.dumps(rsids_grouped_by_ld_prefixes, indent=4))