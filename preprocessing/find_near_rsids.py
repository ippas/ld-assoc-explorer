import sys
import pandas as pd
import requests

url = "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/search/findByChromBpLocationRange"

rsid = sys.argv[1]
chrom = sys.argv[2]
bp = sys.argv[3]

print(rsid)

#5000000
params = {
    "chrom": chrom,
    "bpStart": (int(bp) - 100),
    "bpEnd": (int(bp) + 100),
    "page": 0,
    "size": 499,
    "sort": "bp,asc"
}

response = requests.get(url, params=params)

if response.status_code == 200:
    data = response.json()
    
    for snp in data["_embedded"]["singleNucleotidePolymorphisms"]:
        rsid = snp["rsId"]
        print(rsid)
else:
    print("Błąd podczas pobierania danych:", response.status_code)
