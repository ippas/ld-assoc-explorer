import sys
import json
import pandas
import scipy.sparse as scipy_sparse
import multiprocessing

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

def linked_snps_by_ld_prefix_and_rsid_list(ld_prefix, rsid_list):
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

            linked_snps.append({
                'RS_ID': rs_id,
                'CHR': chr,
                'BP': bp,
                'BASE_SNP': column_name,
                'LD_VALUE': float(value)
            })
    return linked_snps

def process_rsids_group(rsids_group):
    ld_prefix = rsids_group['ld_prefix']
    linked_snps = linked_snps_by_ld_prefix_and_rsid_list(ld_prefix, rsids_group['rs_ids'])


    return {
        'ld_prefix': ld_prefix,
        'linked_snps': linked_snps
    }

def process_worker(input_queue, result_queue):
    while not input_queue.empty():
        input_value = input_queue.get()
        result = process_rsids_group(input_value)
        result_queue.put(result)

def save_progress(found_snps, completed_groups):
    dataframe = pandas.DataFrame(found_snps)
    dataframe.to_csv(results_file_path, index=False)

    with open(completed_rsids_groups_file_path, 'w') as file:
        json.dump(completed_groups, file, indent=4)

def main():
    rsids_grouped_by_ld_prefixes = read_json_file(grouped_rsids_file_path)
    completed_rsids_groups = read_json_file(completed_rsids_groups_file_path)

    rsids_groups_to_process = [
        rsid_group
        for rsid_group in rsids_grouped_by_ld_prefixes
        if rsid_group['ld_prefix'] not in completed_rsids_groups
    ]

    if len(rsids_groups_to_process) == 0:
        return

    found_snps = array_of_dictionaries_from_csv_file(results_file_path)

    input_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()

    for input in rsids_groups_to_process:
        input_queue.put(input)

    num_of_workers = min(3, multiprocessing.cpu_count(), len(rsids_groups_to_process))
    workers = []

    for _ in range(num_of_workers):
        process = multiprocessing.Process(target=process_worker, args=(input_queue, result_queue))
        workers.append(process)
        process.start()

    while any(process.is_alive() for process in workers) or not result_queue.empty():
        result = result_queue.get()
        completed_rsids_groups.append(result['ld_prefix'])
        for found_snp in result['linked_snps']:
            if (any([
                snp['RS_ID'] == found_snp["RS_ID"] and snp['BASE_SNP'] == found_snp['BASE_SNP']
                for snp in found_snps
                ])):
                continue

            found_snps.append(found_snp)
        
        save_progress(found_snps, completed_rsids_groups)

    for process in workers:
        process.join()

    save_progress(found_snps, completed_rsids_groups)

if __name__ == '__main__':
    main()  