# ld-assoc-explorer (Linkage Disequilibrium Associations Explorer)
This project searches for associations between a list of SNPs and phenotypes, taking into account nearby SNPs with high linkage disequilibrium (LD).

## Methods
The program takes a list of SNPs (containing rsID, chromosome, and position in GRCh37 format) and uses UK Biobank LD matrices to find additional SNPs with high LD (r² ≥ 0.8). It then searches for phenotype associations (filtered at p-value ≤ 1e-6) for both the original and LD-expanded SNP lists using:

1. [EBI GWAS Catalog Summary Statistics API](https://www.ebi.ac.uk/gwas/summary-statistics/docs/)
2. [EBI GWAS Catalog top associations](https://www.ebi.ac.uk/gwas/docs/file-downloads)
3. [Open Targets Platform GraphQL API](https://api.platform.opentargets.org/)

**Data Source**:  
UK Biobank Linkage Disequilibrium Matrices were accessed from [https://registry.opendata.aws/ukbb-ld](https://registry.opendata.aws/ukbb-ld) on [DATE].

## Pipeline

### 1. Preprocessing
Groups SNPs by chromosome and position for efficient retrieval from LD matrices (stored as `<CHR>_<BP1>_<BP2>` files).

**Script**: `./preprocessing/group-snp-by-chr-and-bp.py`  
**Input**: `./data/snps-for-analysis.csv` (copy of user-provided SNP list)  
**Output**: `./data/rsids-grouped-by-ld-prefixes.json`

### 2. Analysis

#### Step 1: LD Expansion
Finds SNPs in high LD with the input SNPs.

**Script**: `./analysis/find-linked-snps-by-ld-matrixes.py`  
**Inputs**: 
- `./data/rsids-grouped-by-ld-prefixes.json`
- `./data/snps-for-analysis.csv`  

**Outputs**: 
- `./data/snps-found-via-ld-matrixes.csv` (expanded SNP list)
- `./data/completed-rsids-groups.json` (checkpoint file)

#### Step 2: Association Search
Queries phenotype associations for the expanded SNP list:

| Script | API Used | Output File | Notes |
|--------|----------|-------------|-------|
| `find-associations-using-gwas-catalog.py` | GWAS Catalog API | `associations-found-by-gwas-catalog.csv` | Slower, works for any p-value |
| `find-associations-using-gwas-top-associations.py` | GWAS Catalog top associations | `associations-found-by-gwas-top-associations.csv` | Faster, but only for p ≤ 1e-5 |
| `find-associations-using-open-targets.py` | Open Targets API | `associations-found-by-open-targets.csv` | |

Final results are copied to the `./results` directory.

## Usage

### Option 1: Direct Execution
```bash
./analysis/run-analysis.sh [OPTIONS] SNPS_FILE LD_MATRIXES_DIRECTORY
```

### Option 2: Docker
1. Build the image from `./Dockerfile`
2. Run with required mounts:
```bash
docker run -v /host/input.csv:/input.csv \
           -v /host/results:/usr/src/app/results \
           -v /host/data:/usr/src/app/data \
           -v /host/ld_matrices:/ld_matrices \
           [IMAGE] [OPTIONS] /input.csv /ld_matrices
```

### Arguments

| Option | Description |
|--------|-------------|
| `--ld-value <0-1>` | LD threshold (default: 0.75) |
| `--p-value <0-1>` | p-value threshold (default: 0.05) |
| `--api <API_LIST>` | APIs to use (comma-separated). Options: `gwas-catalog`, `gwas-top-associations`, `open-targets` (default: all) |
| `--max-tries <N>` | Maximum retry attempts for failed commands |
| `-r/--resume` | Resume interrupted analysis |
| `-s/--skip-ld` | Skip LD expansion step |
| `-h/--help` | Show help |

**Important Notes**:
1. `gwas-top-associations` requires downloading the top associations TSV file from [EBI GWAS](https://www.ebi.ac.uk/gwas/docs/file-downloads) and either:
   - Placing it in the root directory, OR
   - Updating the path in `./.config`
2. For p-values > 1e-5, only `gwas-catalog` will work (slower).
