#!/bin/bash

# =============================
# CONFIGURATION
# =============================
BIN_DIR="./bin"            # Folder containing the executables
CO60_DIR="./co60"          # Output directory for "./co60_linfit" results
PEAKS_DIR="./peaks"        # Output directory for all "./Calibration_HistMaker" results
CAL_FILE="/tig/pterodon_data3/S2276/CalibrationFile_Oct172025.cal" #Calibration file name (edit as needed)
# add all anlysistree root file from 60co as you need:
ANALYSIS_FILES_CO60=(/tig/pterodon_data3/S2276/AnalysisTrees/analysis62062_000*  
                     /tig/pterodon_data3/S2276/AnalysisTrees/analysis62062_001*) 
# Format: "source_name file_path(s)"
SOURCE_LIST=(
  "152eu /tig/pterodon_data3/S2276/AnalysisTrees/analysis62060_000* /tig/pterodon_data3/S2276/AnalysisTrees/analysis62060_001*"
  "133ba /tig/pterodon_data3/S2276/AnalysisTrees/analysis62061_00*"
  "60co ${ANALYSIS_FILES_CO60[*]}"
  "56co /tig/pterodon_data3/S2276/AnalysisTrees/analysis62063_00*"
)
# =============================
# STEP 1: Run co60_linfit
# =============================
echo "===================================="
echo "🔹 STEP 1: Co60 linear fit calibration"
echo "===================================="

if [[ -f "${CO60_DIR}/co60_linfit.dat" ]]; then
  echo "✅ co60_linfit.dat already exists, skipping Co60 fit."
else
  echo "Running Co60 linear fit..."
  "$BIN_DIR/co60_linfit" "$CAL_FILE" "${ANALYSIS_FILES_CO60[@]}"
  echo "Done."

  # Create output directory if missing
  mkdir -p "$CO60_DIR"

  # Move results
  [[ -f "co60_linfit.dat" ]] && mv co60_linfit.dat "$CO60_DIR/"
  [[ -f "co60_linfit.root" ]] && mv co60_linfit.root "$CO60_DIR/"

  echo "Results saved to: $CO60_DIR/"
fi

# =============================
# STEP 2: Run Calibration_HistMaker (parallel version)
# =============================
echo
echo "===================================="
echo "🔹 STEP 2: Build calibration histograms (parallel)"
echo "===================================="

mkdir -p "$PEAKS_DIR"
MAX_PARALLEL=4   # Maximum number of jobs running at the same time

job_count=0

# Loop through each source
for entry in "${SOURCE_LIST[@]}"; do
  source_name=$(echo "$entry" | awk '{print $1}')
  file_paths=$(echo "$entry" | cut -d' ' -f2-)

  echo "------------------------------------"
  echo "[INFO] Launching source: $source_name"
  echo "       Files: $file_paths"
  echo "------------------------------------"

  # Expand all file paths into an array
  files=( $file_paths )
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "⚠️  Warning: No matching files found for $source_name"
    continue
  fi

  (
    # ======= Start background job =======
    # Run Calibration_HistMaker for this source
    "$BIN_DIR/Calibration_HistMaker" "$CAL_FILE" "$source_name" "${files[@]}"

    # Move result files to the peaks folder
    [[ -f "peaks_${source_name}.dat" ]] && mv "peaks_${source_name}.dat" "$PEAKS_DIR/"
    [[ -f "peaks_${source_name}.root" ]] && mv "peaks_${source_name}.root" "$PEAKS_DIR/"

    echo "✅ Finished $source_name"
    # ======= End background job =======
  ) &  # Run this block in the background

  ((job_count++))

  # Limit number of concurrent jobs
  if (( job_count >= MAX_PARALLEL )); then
    wait              # Wait for all running jobs to finish
    job_count=0       # Reset job counter
  fi
done

# Wait for any remaining background jobs to complete
wait

echo "✅ All sources finished. All peak files moved to: $PEAKS_DIR/"

# =============================
# STEP 3: Final Calibration
# =============================
echo
echo "===================================="
echo "🔹 STEP 3: Combine calibration sources"
echo "===================================="

# Gather all sources automatically from Peaks folder
SOURCES=()
for file in "$PEAKS_DIR"/peaks_*.dat; do
  [[ -e "$file" ]] || continue
  base=$(basename "$file" .dat)
  src=${base#peaks_}
  SOURCES+=("$src")
done

if [[ ${#SOURCES[@]} -eq 0 ]]; then
  echo "⚠️  No peak data found in $PEAKS_DIR/. Nothing to calibrate."
  exit 1
fi

echo "Found sources: ${SOURCES[*]}"
echo "Running final calibration..."
"$BIN_DIR/Calibration" "${SOURCES[@]}"

# Final outputs stay in current directory
if [[ -f "calibration.root" && -f "cal_pars.dat" ]]; then
  echo "✅ Calibration completed successfully."
  echo "Files generated:"
  echo "  - calibration.root"
  echo "  - cal_pars.dat"
else
  echo "⚠️  Calibration did not produce expected output files."
fi

echo
echo "===================================="
echo "🎯 ALL CALIBRATION STEPS COMPLETE!"
echo "===================================="

