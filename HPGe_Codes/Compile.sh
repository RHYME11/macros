#!/bin/bash

# =============================
# CONFIGURATION
# =============================
SRC_DIR="src"
BIN_DIR="bin"

# Create bin directory if it doesn't exist
mkdir -p "$BIN_DIR"

# Common compilation flags
COMMON_FLAGS="-Wl,--no-as-needed `root-config --cflags --libs --glibs` \
-lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA \
-L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries \
`grsi-config --cflags --all-libs --GRSIData-libs` \
-I$GRSISYS/GRSIData/include"

# =============================
# COMPILE THREE PROGRAMS
# =============================
echo "===================================="
echo "🚀 Compiling co60_linfit.cxx ..."
g++ "$SRC_DIR/co60_linfit.cxx" $COMMON_FLAGS -o "$BIN_DIR/co60_linfit" || { echo "❌ Failed: co60_linfit"; exit 1; }

echo "===================================="
echo "🚀 Compiling Calibration_HistMaker.cxx ..."
g++ "$SRC_DIR/Calibration_HistMaker.cxx" $COMMON_FLAGS -o "$BIN_DIR/Calibration_HistMaker" || { echo "❌ Failed: Calibration_HistMaker"; exit 1; }

echo "===================================="
echo "🚀 Compiling Calibration.cxx ..."
g++ "$SRC_DIR/Calibration.cxx" $COMMON_FLAGS -o "$BIN_DIR/Calibration" || { echo "❌ Failed: Calibration"; exit 1; }

echo "===================================="
echo "✅ All programs compiled successfully!"
echo "Executables saved to: $BIN_DIR/"

