#!/bin/bash

# ============================================
# Clean script for analysis outputs
# Ask before deleting anything
# ============================================

# Output files / patterns
RAW_HIST="raw_hist.root"
FIT_HIST="fit_hist.root"
HIST="hist.root"
CAL_TXT="Calibration.txt"
RES_CHECK="Res_Check.dat"

# Histogram outputs
HIST_FILES="Hist_*.root"

# Peak / intermediate files
TMP_FILES="fit_output.tmp"

# Directories that may contain outputs
DIRS_TO_CLEAN=("co60" "peaks" "HistFiles")

# --------------------------------------------
# Helper function: ask before deleting
# --------------------------------------------
ask_and_remove() {
  local target="$1"

  if compgen -G "$target" > /dev/null; then
    echo
    echo "⚠️  Found: $target"
    read -p "Delete this? [y/N]: " ans
    case "$ans" in
      y|Y)
        rm -rf $target
        echo "🗑️  Deleted: $target"
        ;;
      *)
        echo "⏩ Skipped: $target"
        ;;
    esac
  fi
}

echo "============================================"
echo " CLEAN ANALYSIS OUTPUT FILES"
echo "============================================"

# --------------------------------------------
# Clean single files
# --------------------------------------------
ask_and_remove "$RAW_HIST"
ask_and_remove "$FIT_HIST"
ask_and_remove "$HIST"
ask_and_remove "$CAL_TXT"
ask_and_remove "$RES_CHECK"
ask_and_remove "$TMP_FILES"

# --------------------------------------------
# Clean histogram outputs
# --------------------------------------------
ask_and_remove "$HIST_FILES"

# --------------------------------------------
# Clean output directories
# --------------------------------------------
for d in "${DIRS_TO_CLEAN[@]}"; do
  if [[ -d "$d" ]]; then
    echo
    echo "⚠️  Directory found: $d/"
    read -p "Delete entire directory $d/? [y/N]: " ans
    case "$ans" in
      y|Y)
        rm -rf "$d"
        echo "🗑️  Deleted directory: $d/"
        ;;
      *)
        echo "⏩ Skipped directory: $d/"
        ;;
    esac
  fi
done

echo
echo "============================================"
echo "✅ Clean finished."
echo "============================================"

