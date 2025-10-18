

Mini .c & .cpp codes for ROOT analysis.


**GausFit_orginal.C:** old version for A = 156 coincidence check;
**GausFit.C:** try organizing codes (half-way), for A = 160 analysis.

## AlphaCalibration.c
Auto calibration code for charged partile detector with triple alpha source;<span style="color:red"> GRSISort Required.</span> </br>
**0. Compiling Commond (1st line in the code txt):** `g++ AlphaCalibration.c -Wl,--no-as-needed `root-config --cflags --libs --glibs` -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm -O2 -Wl,--copy-dt-needed-entries -L/opt/local/lib -lX11 -lXpm `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -o bin/Alphacal`;</br>
**1. Input:** (change it on line 525 & 528) </br>
&nbsp;&nbsp;&nbsp;&nbsp;**1.a FragmentTree + CalibrationFile + starting CH + ending CH;**</br>
&nbsp;&nbsp;&nbsp;&nbsp;**1.b AnalysisTree + CalibrationFile + starting CH + ending CH (for TH1 *tdiff);**</br>
**2. Calibration Math Formula: Pu+Am+Cm.** </br>
2.a Edit lin 328 to `TF1 *fc = tasf(hist, Form("fc_CH%i",ich), min,max,"cl");` for **Gd+Th+Cm**;</br>
**3. Output:**</br>
3.a **Calibration.txt**: includes two array, gain and offset;</br>
3.b **Hist.root:** includes calibrated summary TH2 for calibration quick check;</br>
3.c **Values Printed on Screen:** FWHM of three alpha peaks after calibration.</br>
<span style="color:red">Return FWHM = -1, GAIN = 1, OFFSET = 0, **if the channel is empty!**</span>

# Lifetime Measurement Codes
Codes in this section work for lifetime measurement analysis.

## CalChi2.cxx
**Calculate the chi2 for the lineshape shift between experimental data and simulation data. Return chi2 with simulated lifetime-value.** </br>
- Compiling command: line 1 </br>
- Input: lower gamma energy, higher gamma, data root file, simulation root file. </br>
- Notes: </br>
&nbsp;&nbsp;&nbsp;&nbsp; 1. histogram names need to be justed (line68~69); </br>
&nbsp;&nbsp;&nbsp;&nbsp; 2. binwidth for histograms need to be justed (line73~74); </br>
&nbsp;&nbsp;&nbsp;&nbsp; 3. scaling factor: line 79; </br>
&nbsp;&nbsp;&nbsp;&nbsp; 4. all output will print out on the screen; </br>


## runCalChi2.sh
**Run Calchi2.cxx with multiple simulation files with different lifetimes but with the same experimental data file. Will record the minimum chi2 and its corresponding scaling factor for each simulation file. The output is a txt file.** </br>
I didn't upload it but it should be combined to CalChi2.cxx or make it human readable;

## FitChi2.cxx
**The input required the format of runCalChi2.sh output file. Then fit the curve "chi2 vs lifetime". Return min_chi2 and its lifetime. Also get two lifetimes with chi2 = chi2_min+1**  
