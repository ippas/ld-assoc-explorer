#!/bin/bash
#
# This program searches associations for a list of snps taking into account associations with nearby SNPs with high LD. 

scr_dir=$(dirname "$0")
scr_name=$(basename "$0")
cd $scr_dir

# START OF THE SCRIPT

available_apis=('gwas-catalog' 'gwas-top-associations' 'open-targets' 'off')
selected_apis=(${available_apis[@]})

ld_value=0.75
p_value=0.05
snps_file_path=
ld_matrixes_directory=

max_tries=5

do_resume=0
skip_ld=0

analysis_start_time=
last_run_time=

PIDS=()

usage() {
    _usage="\
Usage: $scr_name [OPTION]... SNPS_FILE LD_MATRIXES_DIRECTORY
Searches for associations of provided SNPS from SNPS_FILE, taking into account nearby SNPS with high ld from LD_MATRIXES_DIRECTORY

SNPS_FILE must have .csv extension

  --ld-value <ld_value>               Percentage value between 0 to 1. If not precised: standard value of 0.75
  --p-value <p_value>                 Percentage value between 0 to 1. If not precised: standard value of 0.05
  --api <API_NAME_1,API_NAME_2,...>   Specify one or more API names to connect to, separated by commas (default: all except 'off'). Possible API Names: 
                                        $(echo ${available_apis[@]} | sed 's/ /, /g')
                                      Note: gwas-top-associations works only on p-value<1e-5, for higher p-value use gwas-catalog.
  --max-tries <max_tries>             Maximum number of times the program will try to run a command.
  -r, --resume                        Resumes a previously interrupted analysis using the same data and options as the
                                        prior run, instead of starting a new analysis.
  -s, --skip-ld                       Skips the step of identifying genes with high LD.
                                        This option allows the analysis to proceed without considering finding genes by high LD.
  -h, --help                          Displays this help and exit
"
    echo "$_usage"
    exit 1
}

settings_info() {
    _info="\
====================================
    $scr_name: Info:
====================================
Mode: $( if [ "$do_resume" -eq 1 ]; then echo "Resuming previous analysis"; else echo "Starting new analysis"; fi )\
$( if [ "$do_resume" -eq 1 ]; then echo -e "\nAnalysis Start: $analysis_start_time"; fi )\
$( if [ "$do_resume" -eq 1 ]; then echo -e "\nLast Run: $last_run_time"; fi )

SNPS_FILE           $snps_file_path
$( if [ "$skip_ld" -eq 0 ]; then echo -e "LD_MATRIXES_DIR     $ld_matrixes_directory"; fi )\
$( if [ "$skip_ld" -eq 0 ]; then echo -e "\nLD_VALUE            $ld_value\n \b"; fi )\
P_VALUE             $p_value                                               
API[s]              $(echo ${selected_apis[@]} | sed "s/ /, /g" )          
MAX_TRIES           $max_tries
\
$( if [ "$do_resume" -eq 1 ]; then echo -e "\nNote: Settings cannot be changed in resume mode. Using saved options."; fi )
====================================
"
    echo "$_info"
}

results_info() {
    echo -e "\nAll scripts finished successfully. The output was saved to ../results  "
}

error_message() {
    echo "$scr_name: $1" >&2
}

read_var_state() {
    if [ ! -e "../data/state" ]; then
        error_message "No previous analysis found. A first analysis must be run before resuming."
        exit 1 
    fi

    source "../data/state"
    echo $analysis_start_time
}

save_var_state() {
    _state="\
analysis_start_time=\"$( if [ -e "../data/state" ]; then echo "$analysis_start_time"; else echo "$(date -u +%F\ %T\ UTC)"; fi )\"
last_run_time=\"$(date -u +%F\ %T\ UTC)\"

snps_file_path=$snps_file_path
ld_matrixes_directory=$ld_matrixes_directory

ld_value=$ld_value
p_value=$p_value

selected_apis=(${selected_apis[@]})

max_tries=$max_tries

skip_ld=$skip_ld

"
    echo "$_state" > "../data/state"
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

set_max_tries() {
    if [[ ! "$1" =~ ^[0-9]+$ || "$1" -eq 0 ]]; then
        error_message "--max-tries '$1' has to be positive integer"
        exit 1
    fi

    max_tries=$1
}

kill_all_jobs() {
    for pid in "${PIDS[@]}"; do
        disown "$pid" #2> /dev/null
        kill "$pid" #2> /dev/null
    done
    wait
}

run_command() {
    for try in $(seq $max_tries); do
        command="$@"

        echo "$scr_name: Running command: '$command' ($try/$max_tries)"

        eval "$command"
        exit_code=$?

        if [[ $exit_code -ne 0 && $try -ne $max_tries ]]; then
            error_message "command: '$command' failed and exited with code '$exit_code'"
            continue
        fi

        if [[ $exit_code -ne 0 && $try -eq $max_tries ]]; then
            error_message "command: '$command': maximum tries reached. Exiting program..."
            exit 1
        fi

        echo "$scr_name: command: '$command' finished succesfully"

        return
    done
}

declare_variables() {
    while true; do
        case "$1" in
            --help | -h )
            usage
            ;;
            --resume | -r )
            do_resume=1
            read_var_state
            return
            ;;
            --skip-ld | -s )
            skip_ld=1
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
            --max-tries )
            shift
            set_max_tries $1
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

    if [ "$skip_ld" -eq 0 ]; then
        set_ld_matrixes_directory $2
    fi
}

main() {
    if [ "$#" -eq 0 ]; then
        usage
    fi

    declare_variables "$@"

    settings_info

    if [ "$do_resume" -eq 0 ]; then
        find ../data -regextype egrep -regex '.*\.(csv|json)$' -delete

        if [ -e "../data/state" ]; then
            rm "../data/state"
        fi

        cp $snps_file_path ../data/snps-for-analysis.csv

        if [ "$skip_ld" -eq 0 ]; then
            run_command "python3 ../preprocessing/group-snps-by-chr-and-bp.py ../data/snps-for-analysis.csv $ld_matrixes_directory"
        fi
    fi

    save_var_state

    if [ "$skip_ld" -eq 0 ]; then
        run_command "python3 find-linked-snps-by-ld-matrixes.py $ld_matrixes_directory $ld_value"
        cp "../data/snps-found-via-ld-matrixes.csv" "../results/"
    else
        cp "../data/snps-for-analysis.csv" "../data/snps-found-via-ld-matrixes.csv"
    fi

    for selected_api in ${selected_apis[@]}; do
        case "$selected_api" in
            off )
            break
            ;;
            gwas-catalog )
            run_command "python3 find-associations-using-gwas-catalog.py $p_value" && \
            cp "../data/associations-found-by-gwas-catalog.csv" "../results/" &
            PIDS+=($!)
            ;;
            gwas-top-associations )
            run_command "python3 find-associations-using-gwas-top-associations.py $p_value" && \
            cp "../data/associations-found-by-gwas-top-associations.csv" "../results/" &
            PIDS+=($!)
            ;;
            open-targets )
            run_command "python3 find-associations-using-open-targets.py $p_value" && \
            cp "../data/associations-found-by-open-targets.csv" "../results/" &
            PIDS+=($!)
            ;;
            * )
            error_message "unhandled api '$selected_api'"
            ;;
        esac
    done

    trap 'kill_all_jobs; exit 1' SIGINT

    for pid in "${PIDS[@]}"; do
        wait "$pid" || {
            kill_all_jobs
            exit 1
        }
    done
    
    results_info
    rm "../data/state"
}

main "$@"

# END OF THE SCRIPT

cd - >/dev/null
exit 0