# Project title (institute-name-subject)

## Methods
The preprocessing of raw SNP 

## Preprocessing
Preprocessing of raw data includes grouping snp identifiers by chromosome.

## Analysis
After preprocessing, we used the UK Biobank Linkage Disequilibrium Matrices available from DATE at https://registry.opendata.aws/ukbb-ld.Details. Using this matrixes we find SNP with high LD (> 0.75) for our list of SNP. The next step was to search for associations for our associated SNPs. To do this, I used the summary statistics api from https://www.ebi.ac.uk/gwas/home. I filtered my data for p-value <= 0.05.

*notes: all files included in the repo need to be referenced, either in README or other .md files. The analysis has to be fully reproducible, in principle the repo should contain code + description of how to run it while data and results kept outside*

## About this template
Directories:
- _root_ - README.md, *.Rproj, general configuration files, etc.
- raw - raw data
- preprocessing - scripts
- data - useful data, created by scripts/tools/preprocessing
- analysis - analysis source code
- results - output ready to present
