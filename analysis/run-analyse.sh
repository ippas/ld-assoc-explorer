#!/bin/bash

usage() {
    echo "Usage: $0 <snps_data_file_path> <ld_matrixes_folder_path> [--ld <ld_value> | <ld_value>] [--pvalue <p_value> | -p <p_value>]"
    exit 1
}

if [ "$#" -lt 2 ]; then
    usage
fi

snps_data_file_path="$1"
ld_matrixes_folder_path="$2"
echo "Raw data path: $snps_data_file_path"
echo "LD matrix folder path: $ld_matrices_folder_path"
shift 2

ld_value=75
p_value=0.05


while [[ $# -gt 0 ]]; do
	key="$1"

	case$key in
		--ld|-l)
			ld_value="$2"
			shift 2
			;;
		--pvalue|-p)
			p_value="$2"
			shift 2
			;;
		-*|--*)
			echo "Unkown option $1"
			usage
			;;
		*)
			break
			;;
	esac
done

if [ -n "${ld_value}" ] && ! [[ "${ld_value}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: --ld/-l must be a numeric value."
    usage
fi

if [ -n "${pvalue_value}" ] && ! [[ "${pvalue_value}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: --pvalue/-p must be a numeric value."
    usage
fi

if [ -n "${ld_value}" ]; then
    echo "LD value: $ld_value"
fi

if [ -n "${pvalue_value}" ]; then
    echo "P-value: $pvalue_value"
fi

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
run_script "$script2" "$ld_matrixes_folder_path" "$ld_value"

script3="find-associations.py"
run_script "$script3" "$p_value"

script4="name-associations.py"
run_script "$script3" "$snp_data_file_path" 

script5="format-result.py"
run_script "$script4" 

script6="delete-duplicates.py"
run_script "$script5"

echo "All scripts have been successfully executed."
