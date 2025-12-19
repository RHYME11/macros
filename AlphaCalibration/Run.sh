#!/bin/bash

# ============================================
# Configuration
# ============================================

BIN_DIR="./bins"

RAW_EXE="${BIN_DIR}/RawHistMaker"
FIT_EXE="${BIN_DIR}/FitRawHist"
HIST_EXE="${BIN_DIR}/HistMakers"

CAL_FILE="/data1/yzhu/Projects/S2403/CalibrationFileClean.cal"

# AnalysisTree input (glob will expand automatically)
ANALYSIS_FILES=/data1/yzhu/Projects/S2403/AnalysisTrees/analysis62336*

# Output files
RAW_HIST="raw_hist.root"
RES_CHECK="Res_Check.dat"
TMP_OUT="fit_output.tmp"

# Optional: run HistMakers (1 = yes, 0 = no)
RUN_HISTMAKERS=0

# ============================================
# Step 1: Run RawHistMaker
# ============================================
echo "============================================"
echo "[STEP 1] Running RawHistMaker"
echo "Calibration file : $CAL_FILE"
echo "Analysis files   : $ANALYSIS_FILES"
echo "============================================"

#"$RAW_EXE" "$CAL_FILE" $ANALYSIS_FILES

echo "[OK] Raw histogram created: $RAW_HIST"
echo

# ============================================
# Step 2: Run FitRawHist
# Capture stdout and generate Res_Check.dat
# ============================================
echo "============================================"
echo "[STEP 2] Running FitRawHist"
echo "============================================"

"$FIT_EXE" "$RAW_HIST" > "$TMP_OUT"

awk '
BEGIN {
  OFS = "\t";
}
NR == 1 {
  print "#CHANNEL", "FWHM(Pu)", "FWHM(Am)", "FWHM(Cm)",
        "Res%(5.8MeV)", "GAIN", "OFFSET";
  next;
}
NF >= 6 {
  chan = $1;
  fpu  = $2;
  fam  = $3;
  fcm  = $4;
  gain = $5;
  offs = $6;

  res = fcm / 5800.0 * 100.0;

  printf "%.0f\t%.4f\t%.4f\t%.4f\t%.2f\t%.4f\t%.4f\n",
        chan, fpu, fam, fcm, res, gain, offs;
}
' "$TMP_OUT" > "$RES_CHECK"

rm -f "$TMP_OUT"

echo "[OK] Res_Check.dat created."
echo

# ============================================
# Step 3: Run HistMakers (optional)
# ============================================
if [[ $RUN_HISTMAKERS -eq 1 ]]; then
  echo "============================================"
  echo "[STEP 3] Running HistMakers"
  echo "============================================"

  "$HIST_EXE" "$CAL_FILE" $ANALYSIS_FILES
else
  echo "[INFO] HistMakers step skipped."
fi

echo
echo "============================================"
echo "✅ All steps finished."
echo "Generated files:"
echo "  - $RAW_HIST"
echo "  - Calibration.txt"
echo "  - $RES_CHECK"
echo "============================================"

