#!/bin/bash

run_script() {
    script_name="$1"
    if [ "$script_name" == "find-files-and-find-snps.py" ]; then
        ld_matrix_path="$2"  
        echo "Running script: $script_name with LD matrix path: $ld_matrix_path"
        python "$script_name" "$ld_matrix_path"
    else
        echo "Running script: $script_name"
        python "$script_name"
    fi
    return $?
}

script1="group-rsids-by-chromosome.py"

#in ld_matrices_folder_path are path for folder with LD-matrixes from UK biobank
script2="find-files-and-find-snps.py"
ld_matrixes_folder_path = "/mnt/droplet"

script3="name-associations.py"

script4="format-result.py"

script5="delete-duplicates.py"

scripts=("$script1" "$script2" "$script3" "$script4" "$script5")

for script in "${scripts[@]}"; do
    while true; do
        if [ "$script" == "$script2" ]; then
            run_script "$script" "$ld_matrixes_folder_path"
        else
            run_script "$script"
        fi
        exit_code=$?
        if [ $exit_code -eq 0 ]; then
            break
        else
            echo "Error: Script $script did not finish successfully. Retrying..."
        fi
    done
done

echo "All scripts have been successfully executed."
