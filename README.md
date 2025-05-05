

Mini .c & .cpp codes for ROOT analysis.


**GausFit_orginal.C:** old version for A = 156 coincidence check;
**GausFit.C:** try organizing codes (half-way), for A = 160 analysis.

## AlphaCalibration.c
Auto calibration code for charged partile detector with triple alpha source;<span style="color:red"> GRSISort Required.</span>
**0. Compiling Commond (1st line in the code txt):** `g++ AlphaCalibration.c -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -o bin/Alphacal`;
**1. Input: FragmentTree + CalibrationFile + starting CH + ending CH;**
**2. Calibration Math Formula: Pu+Am+Cm.** 
2.a Edit lin 328 to `TF1 *fc = tasf(hist, Form("fc_CH%i",ich), min,max,"cl");` for **Gd+Th+Cm**;
**3. Output:**
3.a **Calibration.txt**: includes two array, gain and offset;
3.b **Hist.root:** includes calibrated summary TH2 for calibration quick check;
3.c **Values Printed on Screen:** FWHM of three alpha peaks after calibration.
<span style="color:red">Return FWHM = -1, GAIN = 1, OFFSET = 0, **if the channel is empty!**</span>
