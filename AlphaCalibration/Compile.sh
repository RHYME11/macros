#!/bin/bash

# =============================
# Build script for 3 programs
# =============================

CXX=g++
CXXFLAGS="-O2 -Wl,--no-as-needed -Wl,--copy-dt-needed-entries"

# --- Output directory ---
BINDIR="bins"
mkdir -p "$BINDIR"

# --- ROOT flags ---
ROOTFLAGS=$(root-config --cflags --libs --glibs)

# --- GRSI flags ---
GRSIFLAGS=$(grsi-config --cflags --all-libs --GRSIData-libs)

# --- Extra libs ---
EXTRA_LIBS="-lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -lX11 -lXpm -lROOTTPython"

# --- Include paths ---
INCLUDES="-I$GRSISYS/GRSIData/include"

# =============================
# Compile RawHistMaker
# =============================
echo "Compiling RawHistMaker..."
$CXX src/RawHistMaker.cxx $CXXFLAGS \
    $ROOTFLAGS $EXTRA_LIBS \
    -L/opt/local/lib -lX11 -lXpm \
    $GRSIFLAGS $INCLUDES \
    -o "$BINDIR/RawHistMaker"

echo "✔ RawHistMaker built → $BINDIR/RawHistMaker"

# =============================
# Compile FitRawHist
# =============================
echo "Compiling FitRawHist..."
$CXX src/FitRawHist.cxx $CXXFLAGS \
    $ROOTFLAGS $EXTRA_LIBS \
    -L/opt/local/lib -lX11 -lXpm \
    $GRSIFLAGS $INCLUDES \
    -o "$BINDIR/FitRawHist"

echo "✔ FitRawHist built → $BINDIR/FitRawHist"

# =============================
# Compile HistMakers
# =============================
echo "Compiling HistMakers..."
$CXX src/HistMakers.cxx $CXXFLAGS \
    $ROOTFLAGS $EXTRA_LIBS \
    -L/opt/local/lib -lX11 -lXpm \
    $GRSIFLAGS $INCLUDES \
    -o "$BINDIR/HistMakers"

echo "✔ HistMakers built → $BINDIR/HistMakers"

echo "=============================="
echo "✔ All programs compiled into $BINDIR/"
echo "=============================="

