import gzip
import shutil
import requests
import urllib.request
import subprocess

def find_ld_matrix(reference_snps, plink_folder, bfile_prefix):
    plink_ld = [
        "plink",
        "--bfile", f"{plink_folder}/{bfile_prefix}",
        "--window", 0, 999999999,
    ]
    for snp in reference_snps:
        plink_ld.extend(["--ld-snp", snp])
    plink_ld.append("--r2")

    try:
        subprocess.run(plink_ld, check=True)
    except subprocess.CalledProcessError as e:
        print("Błąd podczas wykonania polecenia PLINK:", e)

#vcf_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
#tbi_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"

vcf_file = "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
tbi_file = "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"

#urllib.request.urlretrieve(vcf_url, vcf_file)
#urllib.request.urlretrieve(tbi_url, tbi_file)
plink_folder = "../data"
bfile_prefix = "XD"

plink_vcf_cmd = f"gzip -dc {vcf_file} | plink --vcf /dev/stdin --make-bed --out {plink_folder}/{bfile_prefix}"


try:
    subprocess.run(plink_vcf_cmd, check=True)
except subprocess.CalledProcessError as e:
    print("Błąd podczas wykonania polecenia PLINK:", e)

snps = ["rs189257163", "rs114541806"]
find_ld_matrix(snps, plink_folder, bfile_prefix)
