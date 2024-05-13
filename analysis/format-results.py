import csv
import os
import re
import requests

result_file_path = '../results/result.csv'
associations_path = '../data/snp/associations'
matrixes_path = f'../data/snp/highest-ld-75'
pattern = re.compile(r'rs\d+\.csv')

with open(result_file_path, 'w', newline='') as result:
    writer = csv.DictWriter(result, fieldnames=['chr','main_rsid','ld','second_rsid','efo','trait_name','p_value','study'])
    for root, dirs, files in os.walk(associations_path):
        for file in files:
            if pattern.match(file):
                main_rsid = file.split('.')[0]
                file_path = os.path.join(associations_path, file)
                with open(file_path, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        chrom = row['chromosome']
                        second_rsid = row['rsid']
                        efo = row['trait'][2:-2]
                        study = row['study_accession']
                        p_value = row['p_value']
                        trait_name = efo
                        
                        url= f'https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{efo}'
                        response = requests.get(url)
                        if response.status_code == 200:
                            response = response.json()
                            trait_name = response['trait']
                        
                        matrix_data_path = os.path.join(matrixes_path, f'chrom-{chrom}', f'{main_rsid}.csv')
                        with open(matrix_data_path, 'r') as matdata:
                            matdata_reader = csv.DictReader(matdata)
                            print(matrix_data_path)
                            for matdata_row in matdata_reader:
                                if matdata_row['SNP'] == second_rsid:
                                    bp_allel = f'{chrom}.{matdata_row["BP"]}.{matdata_row["A1"]}.{matdata_row["A2"]}'
                                    matrix_path = os.path.join(matrixes_path, f'chrom-{chrom}', f'{main_rsid}-matrix.csv')
                                    with open(matrix_path, 'r') as mat:
                                        mat_reader = csv.DictReader(mat)
                                        for mat_row in mat_reader:
                                            for col in mat_reader.fieldnames:
                                                if col == bp_allel and col != '':
                                                    ld_value = mat_row[bp_allel]
                                                    break
                                        break
                            writer.writerow({'chr': chrom, 'main_rsid': main_rsid, 'ld': ld_value, 'second_rsid': second_rsid, 'efo': efo, 'trait_name': trait_name, 'p_value': p_value, 'study': study})
