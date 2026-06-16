#!/usr/bin/env bash
#
# run_plots.sh
# Wrapper script to execute cosmological simulation plotting utilities with synchronized arguments.
# Supports isolated inputs for 3D binary snapshots, 1D evolutionary text logs, and metallicity toggles.

# Exit immediately if a command exits with a non-zero status.
set -e

# Define utility directory path
UTILS="/scratch/WORK/SPACE/WORK/Hands_On/Utilities"

# Enforce output workspace presence
mkdir -p plots

# --- Default Configurations ---
LABEL=""
TITLE=""
VOLUME="125.0"
SFR_FILE="sfr.txt"  # Explicit default for the star formation history log
METALS_PRESENT="1"  # Default switch value: 1 = perform metal analysis, 0 = skip

# Phase space default limits (rho and T)
XMIN_RHO="1e-11"
XMAX_RHO="5e-2"
YMIN_TEMP="5e2"
YMAX_TEMP="2e7"

# Star formation history default limits (z and SFR)
XMIN_SFR="0"
XMAX_SFR="12"
YMIN_SFR="0.0001"
YMAX_SFR="0.5"

# Usage/Help signature
usage() {
    cat << EOF
Usage: $0 [options] <snapshot_file_1> [snapshot_file_2]

Options:
  -f FILE   Path to the Star Formation History log file (default: sfr.txt)
  -Z INT    Toggle metallicity analysis. 1=Run metal plots, 0=Skip them (default: 1)
  -l STR    Label for the dataset/legend (default: empty)
  -t STR    Custom title for the plots (default: empty)
  -v NUM    Cosmological simulation volume (default: 125.0)
  -r MIN MAX Density (rho) x-limits (default: 1e-11 1e-4)
  -T MIN MAX Temperature (T) y-limits (default: 5e2 1e7)
  -z MIN MAX Redshift (z) x-limits for SFR history (default: 0 6)
  -s MIN MAX SFR y-limits for SFR history (default: 0 0.8)
  -h        Display this help message
EOF
    exit 1
}

# --- Parse Arguments ---
while getopts "f:Z:l:t:v:r:T:z:s:h" opt; do
    case "${opt}" in
        f) SFR_FILE="${OPTARG}" ;;
        Z) METALS_PRESENT="${OPTARG}" ;;
        l) LABEL="${OPTARG}" ;;
        t) TITLE="${OPTARG}" ;;
        v) VOLUME="${OPTARG}" ;;
        r)
            XMIN_RHO="${OPTARG}"
            eval next_arg=\${$OPTIND}
            if [[ -n "$next_arg" && "$next_arg" != -* ]]; then
                XMAX_RHO="$next_arg"
                OPTIND=$((OPTIND + 1))
            fi
            ;;
        T)
            YMIN_TEMP="${OPTARG}"
            eval next_arg=\${$OPTIND}
            if [[ -n "$next_arg" && "$next_arg" != -* ]]; then
                YMAX_TEMP="$next_arg"
                OPTIND=$((OPTIND + 1))
            fi
            ;;
        z)
            XMIN_SFR="${OPTARG}"
            eval next_arg=\${$OPTIND}
            if [[ -n "$next_arg" && "$next_arg" != -* ]]; then
                XMAX_SFR="$next_arg"
                OPTIND=$((OPTIND + 1))
            fi
            ;;
        s)
            YMIN_SFR="${OPTARG}"
            eval next_arg=\${$OPTIND}
            if [[ -n "$next_arg" && "$next_arg" != -* ]]; then
                YMAX_SFR="$next_arg"
                OPTIND=$((OPTIND + 1))
            fi
            ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND -1))

# --- Input Validation ---
if [ $# -eq 0 ]; then
    echo "Error: Missing mandatory snapshot target file(s)." >&2
    usage
fi

if [ ! -f "$SFR_FILE" ]; then
    echo "Error: Star formation history file '$SFR_FILE' not found." >&2
    exit 1
fi

SNAPSHOT_FILES=("$@")

# Ensure output goes to the created folder if not overridden by structural routines
COMMON_OPTS=()
if [ -n "$LABEL" ]; then
    COMMON_OPTS+=( "--labels" "$LABEL" )
fi
if [ -n "$TITLE" ]; then
    COMMON_OPTS+=( "--title" "$TITLE" )
fi

echo "========================================================================"
echo " Starting Data Visualization Pipeline"
echo "========================================================================"
echo "Snapshots:       ${SNAPSHOT_FILES[*]}"
echo "SFR Log File:    $SFR_FILE"
echo "Metallicity (Z): $METALS_PRESENT"
echo "Volume:          $VOLUME"
echo "Rho Limits:      $XMIN_RHO to $XMAX_RHO"
echo "Temp Limits:     $YMIN_TEMP to $YMAX_TEMP"
echo "SFR z-Limits:    $XMIN_SFR to $XMAX_SFR"
echo "SFR Limits:      $YMIN_SFR to $YMAX_SFR"
echo "------------------------------------------------------------------------"

# --- Command Execution Pipeline ---

# 1. Plot Phase Space Diagram (Density vs Temperature)
#    We modify the internal scaling array colorcode to track mass if metals are skipped.
echo "--> Generating Density-Temperature Phase Space Maps..."
PHASE_OPTS=("${COMMON_OPTS[@]}")
if [ "$METALS_PRESENT" = "0" ]; then
    PHASE_OPTS+=( "--colorcode" "mass" )
fi

echo "--> Generating Phase Space Maps..."
python3 "${UTILS}/plot_phasespace_with_massdistribution.py" \
	"${SNAPSHOT_FILES[@]}" \
	--colorcode=mass \
	--xlim "$XMIN_RHO" "$XMAX_RHO" \
	--ylim "$YMIN_TEMP" "$YMAX_TEMP" \
	"${PHASE_OPTS[@]}"

# 2. Plot Phase Space Diagram with 1D Marginal Mass Distributions
echo "--> Generating Joint Marginals Phase Space Maps..."
python3 "${UTILS}/plot_phasespace_with_massdistribution.py" \
	"${SNAPSHOT_FILES[@]}" \
	--colorcode=mass \
	--xlim "$XMIN_RHO" "$XMAX_RHO" \
	--ylim "$YMIN_TEMP" "$YMAX_TEMP" \
	--marginals \
	"${PHASE_OPTS[@]}"

# Conditional Blocks for Metal-Dependent Plots
if [ "$METALS_PRESENT" = "1" ]; then

    # 3. Plot Phase Space Diagram 
    echo "--> Generating Phase Space Maps for metals..."
    python3 "${UTILS}/plot_phasespace_with_massdistribution.py" \
	    "${SNAPSHOT_FILES[@]}" \
	    --colorcode=metals \
	    --xlim "$XMIN_RHO" "$XMAX_RHO" \
	    --ylim "$YMIN_TEMP" "$YMAX_TEMP" \
	    "${PHASE_OPTS[@]}"
    
    # 4. Plot Phase Space Diagram 
    echo "--> Generating Joint Marginals Phase Space Maps for metals..."
    python3 "${UTILS}/plot_phasespace_with_massdistribution.py" \
	    "${SNAPSHOT_FILES[@]}" \
	    --colorcode=metals \
	    --xlim "$XMIN_RHO" "$XMAX_RHO" \
	    --ylim "$YMIN_TEMP" "$YMAX_TEMP" \
	    --marginals \
	    "${PHASE_OPTS[@]}"

    # 5. Plot Phase Space Diagram 
    echo "--> Generating Phase Space Maps for metallicity..."
    python3 "${UTILS}/plot_phasespace_with_massdistribution.py" \
	    "${SNAPSHOT_FILES[@]}" \
	    --colorcode=metallicity \
	    --xlim "$XMIN_RHO" "$XMAX_RHO" \
	    --ylim "$YMIN_TEMP" "$YMAX_TEMP" \
	    "${PHASE_OPTS[@]}"
    
    # 6. Plot Phase Space Diagram 
    echo "--> Generating Joint Marginals Phase Space Maps for metallicity..."
    python3 "${UTILS}/plot_phasespace_with_massdistribution.py" \
	    "${SNAPSHOT_FILES[@]}" \
	    --colorcode=metallicity \
	    --xlim "$XMIN_RHO" "$XMAX_RHO" \
	    --ylim "$YMIN_TEMP" "$YMAX_TEMP" \
	    --marginals \
	    "${PHASE_OPTS[@]}"
    
    # 7. Plot metallicity distribution for gas and stars
    echo "--> Generating stars and gas Metallicity Distribution.."
    python3 "${UTILS}/plot_metallicity_distribution.py" \
            "${SNAPSHOT_FILES[@]}" \
            -q Z_both \
            "${COMMON_OPTS[@]}"
    
    # 8. Plot metallicity distribution for gas by temperature
    echo "--> Generating gas Metallicity Distribution by Temperature.."
    python3 "${UTILS}/plot_metallicity_distribution.py" \
            "${SNAPSHOT_FILES[@]}" \
            -q T_gas_metals \
            "${COMMON_OPTS[@]}"
    
    # 9. Plot metallicity distribution for gas by Density
    echo "--> Generating gas Metallicity Distribution by Density.."
    python3 "${UTILS}/plot_metallicity_distribution.py" \
            "${SNAPSHOT_FILES[@]}" \
            -q rho_gas_metals \
            "${COMMON_OPTS[@]}"
else
    echo "--> Skip-Flag [-Z 0] Detected: Bypassing 1D Metallicity Distribution Plots."
fi

# 10. Plot Star Formation History
echo "--> Generating Star Formation Rate History Plots..."
python3 "${UTILS}/plot_sfr.py" \
	"$SFR_FILE" \
	-d \
	-V "$VOLUME" \
	--xlim "$XMIN_SFR" "$XMAX_SFR" \
	--ylim1 "$YMIN_SFR" "$YMAX_SFR" \
	"${COMMON_OPTS[@]}"

echo "========================================================================"
echo " Run Completed Successfully."
echo "========================================================================"

mv -f *png plots
