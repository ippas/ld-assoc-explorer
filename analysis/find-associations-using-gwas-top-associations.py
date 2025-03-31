import sys
import pandas
import math

P_VALUE = float(sys.argv[1])
P_VALUE_MLOG = math.log10(P_VALUE) * -1
SNP_TO_ANALISE = '../data/snps-found-via-ld-matrixes.csv'
RESULTS_FILE_PATH = '../data/associations-found-by-gwas-top-associations.csv'
CONFIG_FILE='../.config'

def read_config_file():
    config = {}

    with open(CONFIG_FILE, 'r') as f:
        for line in f:
            line_formated = line.rstrip()

            if line_formated == '':
                continue

            if line_formated[0] == '#':
                continue

            config[line_formated.split('=')[0]] = line_formated.split('=')[1]
    return config

def main():
    snps_to_analise = pandas.read_csv(SNP_TO_ANALISE)

    config = read_config_file()

    associations_file = pandas.read_csv(config['gwas_catalog_association_file_path'], sep='\t', low_memory=False)

    associations_filtered = associations_file.loc[
                                    associations_file['PVALUE_MLOG'] >= P_VALUE_MLOG, 
                                    ['SNPS', 'PVALUE_MLOG', 'OR or BETA', 'PUBMEDID', 'STUDY ACCESSION', 'MAPPED_TRAIT_URI', 'DISEASE/TRAIT']
                                ]

    associations_filtered = associations_filtered[associations_file['SNPS'].isin(snps_to_analise['RS_ID'])]

    merged_associations_and_snps = snps_to_analise.merge(associations_filtered, left_on='RS_ID', right_on='SNPS', how='left')

    merged_associations_and_snps = merged_associations_and_snps.rename(columns={
        'STUDY ACCESSION': 'STUDY_ID',
        'OR or BETA': 'OR/BETA'
        })

    merged_associations_and_snps.to_csv(RESULTS_FILE_PATH, index=False)

if __name__ == '__main__':
    main()