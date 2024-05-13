#!/bin/bash

#in snps_data_file_path are filepath to our raw dataset
snps_data_file_path="../raw/nearest_with_rare_march.csv"
echo "Raw data path: $snps_data_file_path"

#in ld_matrixes_folder_path are path for folder with LD-matrixes from UK biobank
ld_matrixes_folder_path="/mnt/droplet"
echo "LD matrix folder path: $ld_matrices_folder_path"

run_script() {
    script=$1
    shift 1
    echo "Running script: $script"
    sleep 5
    success=false
    while ! $success; do 
        python "$script" "$@"
        exit_code=$?
        if [ $exit_code -eq 0 ]; then
            success=true
        else
            echo "Error: Script $script did not finish successfully. Retrying..."
        fi
    done
}

script1="../preprocessing/group-rsids-by-chromosome.py"
run_script "$script1" "$snps_data_file_path"

script2="find-files-and-find-snps.py"
run_script "$script2" "$ld_matrixes_folder_path"

script3="name-associations.py"
run_script "$script3" 

script4="format-result.py"
run_script "$script4" 

script5="delete-duplicates.py"
run_script "$script5"

echo "All scripts have been successfully executed."
