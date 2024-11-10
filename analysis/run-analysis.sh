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
available_apis=('gwas-catalog' 'iue-gwas' 'open-targets')
selected_apis=(${available_apis[@]})
is_quiet=0

usage() {
    echo "Usage: $scr_name [OPTION]... SNPS_FILE LD_MATRIXES_DIRECTORY                                                                                          "
    echo "Searches for associations of provided SNPS from SNPS_FILE, taking into account nearby SNPS with high ld from LD_MATRIXES DIRECTORY                    "
    echo "SNPS_FILE must have .csv extension                                                                                                                    "
    echo "                                                                                                                                                      "
    echo "  --help | -h                               Displays this help and exit                                                                               "
    echo "  --quiet | -q                              Runs program in quiet mode - hides unneccesary information.                                               "
    echo "  --ld-value <ld_value>                     Percentage value between 0 to 1. If not precised: standard value of 0.75                                  "
    echo "  --p-value <p_value>                       Percentage value between 0 to 1. If not precised: standard value of 0.05                                  "
    echo "  --api <API_NAME_1,API_NAME_2,...>         Specify one or more API names to connect to, separated by commas (default: all). Possible API Names:      "
    echo "                                              $(echo ${available_apis[@]} | sed 's/ /, /g')                                                           "
    exit 1
}

settings_info() {
    echo "$scr_name: Settings:"
    echo "  SNPS_FILE           $snps_file_path                                        "
    echo "  LD_MATRIXES_DIR     $ld_matrixes_directory                                 "
    echo "  LD_VALUE            $ld_value                                              "
    echo "  P_VALUE             $p_value                                               "
    echo "  API[s]              $(echo ${selected_apis[@]} | sed "s/ /, /g" )          "
    echo "                                                                             "
}

results_info() {
    echo -e "\n$scr_name: All scripts finished successfully. The output was saved to ../results  "
}

error_message() {
    echo "$scr_name: $1" >&2
}

set_ld_value() {
    if [[ ! "$1" =~ ^((0(\.[0-9]+)?)|1)$ ]]; then
        error_message "--ld-value '$1' has to be float in range 0 to 1"
        exit 1
    fi

    ld_value=$1
}

set_p_value() {
    if [[ ! "$1" =~ ^((0(\.[0-9]+)?)|1)$ ]]; then
        error_message "--p-value '$1' has to be float in range 0 to 1"
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

set_apis() {
    selected_apis=()
    for api_arg in $(echo "$1" | sed 's/,/ /g' ); do
        if [[ " ${selected_apis[@]} " =~ " $api_arg " ]]; then
            continue
        fi
        if [[ ! " ${available_apis[@]} " =~ " $api_arg " ]]; then
            error_message "unknown API name '$api_arg'"
            exit 1
        fi

        selected_apis+=($api_arg)
    done
}

run_command() {
    max_tries=1

    for try in $(seq $max_tries); do
        command="$@"

        if [ "$is_quiet" -eq 0 ]; then
            echo "$scr_name: Running command: '$command' ($try/$max_tries)"
        fi

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
        --quiet | -q )
        is_quiet=1
        ;;
        --ld-value )
        shift
        set_ld_value $1
        ;;
        --p-value ) 
        shift
        set_p_value $1
        ;;
        --api )
        shift
        set_apis $1
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

if [ "$is_quiet" -eq 0 ]; then
    settings_info
fi

run_command find ../data -regextype egrep -regex '.*(\.csv|\.json)' -delete

run_command python3 ../preprocessing/group-snps-by-chr-and-bp.py $snps_file_path $ld_matrixes_directory
run_command python3 find-linked-snps-by-ld-matrixes.py $ld_matrixes_directory $ld_value

for selected_api in ${selected_apis[@]}; do
    case "$selected_api" in
        gwas-catalog )
        run_command python3 find-associations-using-gwas-catalog.py $p_value
        ;;
        iue-gwas )
        echo 'iue-gwas'
        ;;
        open-targets )
        run_command python3 find-associations-using-open-targets.py
        ;;
        * )
        error_message "unhandled api '$selected_api'"
        ;;
    esac
done

if [ "$is_quiet" -eq 0 ]; then
    results_info
fi

# END OF THE SCRIPT

cd - >/dev/null
exit 0