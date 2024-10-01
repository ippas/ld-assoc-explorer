#!/bin/bash
#
# This program searches associations for a list of snps taking into account associations with nearby SNPs with high LD. 

scr_dir=$(dirname "$0")
scr_name=$(basename "$0")
cd $scr_dir

# START OF THE SCRIPT

ld_value=0.75
p_value=0.05
snps_file_path=
ld_matrixes_directory=

usage() {
    echo "Usage: $scr_name [OPTION]... SNPS_FILE LD_MATRIXES_DIRECTORY                                                                           "
    echo "Searches for associations of provided SNPS from SNPS_FILE, taking into account nearby SNPS with high ld from LD_MATRIXES DIRECTORY     "
    echo "SNPS_FILE must have .csv extension                                                                                                     "
    echo "                                                                                                                                       "
    echo "  --ld-value <ld_value>       Percentage value between 0 to 1. If not precised: standard value of 0.75                                 "
    echo "  --p-value <p_value>         Percentage value between 0 to 1. If not precised: standard value of 0.05                                 "
    echo "  --help | -h                 Displays this help and exit                                                                              "
    exit 1
}

settings_info() {
    echo "$scr_name: Settings:"
    echo "  SNPS_FILE           $snps_file_path         "
    echo "  LD_MATRIXES_DIR     $ld_matrixes_directory  "
    echo "  LD_VALUE            $ld_value               "
    echo "  P_VALUE             $p_value                "
    echo "                                              "
}

results_info() {
    echo -e "\n$scr_name: All scripts finished successfully. The output was saved to ../results  "
}

error_message() {
    echo "$scr_name: $1" >&2
}

set_ld_value() {
    if [[ ! "$1" =~ ^((0(\.[0-9]+)?)|1)$ ]]; then
        echo "$scr_name: --ld-value '$1' has to be float in range 0 to 1"
        exit 1
    fi

    ld_value=$1
}

set_p_value() {
    if [[ ! "$1" =~ ^((0(\.[0-9]+)?)|1)$ ]]; then
        error_message "$--p-value '$1' has to be float in range 0 to 1"
        exit 1
    fi

    p_value=$1
}

set_snps_file_path() {
    if [ -z "$1" ]; then
        error_message "Missing SNPS file name"
        exit 1
    fi

    if [[ ! (-f "$1" && "$1" =~ .csv$) ]]; then
        error_message "$1: No such SNPS file or it isn't in .csv extension"
        exit 1
    fi

    snps_file_path="$1"
}

set_ld_matrixes_directory() {
    if [ -z "$1" ]; then
        error_message "Missing ld matrixes directory name"
        exit 1
    fi

    if [[ ! -d "$1" ]]; then
        error_message "$1: The provided ld matrixes directory doesn't exist"
        exit 1
    fi

    ld_matrixes_directory="$1"
}

run_command() {
    max_tries=1

    for try in $(seq $max_tries); do
        command="$@"

        echo "$scr_name: Running command: '$command' ($try/$max_tries)"

        $command
        exit_code=$?

        if [[ $exit_code -ne 0 && $try -ne $max_tries ]]; then
            error_message "command: '$command' failed and exited with code '$exit_code'"
            continue
        fi

        if [[ $exit_code -ne 0 && $try -eq $max_tries ]]; then
            error_message "command: '$command': maximum tries reached. Exiting program..."
            exit 1
        fi

        return
    done
}


if [ "$#" -eq 0 ]; then
    usage
fi

while true; do
    case "$1" in
        --help | -h )
        usage
        ;;
        --ld-value )
        shift
        set_ld_value $1
        ;;
        --p-value ) 
        shift
        set_p_value $1
        ;;
        -* )
        error_message "invalid option '$1'"
        exit 1
        ;;
        * )
        break
        ;;
    esac
    shift
done

set_snps_file_path $1
set_ld_matrixes_directory $2

settings_info

run_command python3 ../preprocessing/snps-to-ld-matrixes-format.py $snps_file_path $ld_matrixes_directory "../data/snps-ld-matrixes-format.json"

results_info

# END OF THE SCRIPT

cd - >/dev/null
exit 0