#!/bin/bash
# =============================================
# Clean-up Script for Calibration Environment
# =============================================

CO60_DIR="./co60"
PEAKS_DIR="./peaks"

echo "==============================================="
echo "🧹 Cleaning calibration workspace"
echo "==============================================="

# 1️⃣  Remove all *.root and *.dat files in the current directory
echo "[INFO] Removing all *.root and *.dat files in current directory..."
rm -f ./*.root ./*.dat
echo "✅ Done."

# 2️⃣  Clean peaks/ folder (no confirmation)
if [[ -d "$PEAKS_DIR" ]]; then
  echo "[INFO] Cleaning $PEAKS_DIR/ ..."
  rm -rf "$PEAKS_DIR"/*
  echo "✅ Peaks directory cleaned."
else
  echo "[INFO] Peaks directory not found. Skipping."
fi

# 3️⃣  Ask before cleaning co60/ folder
if [[ -d "$CO60_DIR" ]]; then
  echo
  read -p "⚠️  Do you want to delete all files in $CO60_DIR/? (y/N): " confirm
  if [[ "$confirm" == "y" || "$confirm" == "Y" ]]; then
    rm -rf "$CO60_DIR"/*
    echo "✅ Co60 directory cleaned."
  else
    echo "❎ Skipped cleaning $CO60_DIR/."
  fi
else
  echo "[INFO] Co60 directory not found. Skipping."
fi

echo
echo "🎯 Clean-up completed."

