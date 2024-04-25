#!/bin/bash

main_directory="../data/matrixes/"

total_csv_files=$(find "$main_directory" -name "*.csv")

total_count=$(echo "$total_csv_files" | wc -l)

echo "$total_count"

echo "$total_csv_files" | sed 's/.*\///; s/\.csv$//' > ../data/ieugwasr_list.csv
