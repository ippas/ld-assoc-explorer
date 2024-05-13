# Project title (ifpan-paula-ukbb-ld-matrixes)

This project searches associations for a list of snps taking into account associations with nearby SNPs with high LD.  

## Preprocessing
Preprocessing of raw data includes grouping snp identifiers by chromosome.

## Methods
This sections should be a description of preprocessin and analysis ready to be included in the publication

## Analysis
After preprocessing, we used the UK Biobank Linkage Disequilibrium Matrices available from DATE at https://registry.opendata.aws/ukbb-ld.Details. Using this matrixes we find SNP with high LD (> 0.75) for our list of SNP. The next step was to search for associations for list of associated SNPs. For this purpose, the summary statistics API from [EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/home) was utilized. Data was filtered for a p-value <= 0.05.
To repeat analysis you need your list of snps of interest and folder with LD matrixes from UK BioBank (you need to update filepaths to this data in run-analyse.sh file). To start the analysis, run the run-analyse.sh file.

## About this template
Directories:
- _root_ - README.md, *.Rproj, general configuration files, etc.
- raw - raw data
- preprocessing - scripts
- data - useful data, created by scripts/tools/preprocessing
- analysis - analysis source code
- results - output ready to present
