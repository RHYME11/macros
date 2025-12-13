

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
&nbsp;&nbsp;&nbsp;&nbsp; 1. histogram names need to be justed (line68\~69); </br>
&nbsp;&nbsp;&nbsp;&nbsp; 2. binwidth for histograms need to be justed (line73\~74); </br>
&nbsp;&nbsp;&nbsp;&nbsp; 3. scaling factor: line 79; </br>
&nbsp;&nbsp;&nbsp;&nbsp; 4. all output will print out on the screen; </br>


## runCalChi2.sh
**Run Calchi2.cxx with multiple simulation files with different lifetimes but with the same experimental data file. Will record the minimum chi2 and its corresponding scaling factor for each simulation file. The output is a txt file.** </br>
I didn't upload it but it should be combined to CalChi2.cxx or make it human readable;

## FitChi2.cxx
**The input required the format of runCalChi2.sh output file. Then fit the curve "chi2 vs lifetime". Return min_chi2 and its lifetime. Also get two lifetimes with chi2 = chi2_min+1** 


# HPGe_Codes
The current version can reach based on TIGRESS dataset with AnalysisTree:</br>
1. Resolution check with 60Co; </br>
2. Quadic calibration with different sources for **each crystal**; </br>

## Complie
1. Run `bash Compile.sh` to compile all three codes together; </br>
2. The compiling commands of individual .cxx file are in its 1st line; </br>

## FWHM Check
If you only want to check FWHM, you only need "co60_linfit.cxx": </br>
1. Compile it first; </br>
2. Run `./co60_linfit calibrationfile analysistree_files`; </br>
3. FWHM of each crystal will be saved in "co60_linfit.dat";

## Calibration
1. Edit Run.sh:
&nbsp;&nbsp;&nbsp;&nbsp; 1.1 line9: edit calibration file path; </br>
&nbsp;&nbsp;&nbsp;&nbsp; 1.2 line11\~12: edit co60 analysistree root file paths; </br>
&nbsp;&nbsp;&nbsp;&nbsp; 1.3 line14\~19: edit sources and their relative analysistree root file paths; </br>
2. Run `bash Run.sh`

| Step in Run.sh | .cxx file                 | Input                                                                                  | Output                                                                                                                                                                                                                                                                                                | Notes                                                                                                                                                                                                                                                                                              |
|----------------|---------------------------|----------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Step1          | co60_linfit.cxx           | 1. Calibration_File </br> 2. AnalysisTree root files (you can put multiple root files) | 1. co60_linfit.dat, **including FWHM info** </br> 2. c060_linfit.root:</br>      2.1 uncalibrated histogram of each crystal, the peakfitting results can be reached by: Eg, 1st peak fitting: `TF1 *fx = (TF1 *)hs16->GetListOfFunctions()->At(1)`</br>       2.2 linear calibration of each crystal; | 1. Fit Co60 spectrum with linear function</br> 2. Input must be 60Co relative analysis root files </br>                                                                                                                                                                                            |
| Step2          | Calibration_HistMaker.cxx | 1. Calibration_File </br> 2. Source name </br> 3. AnalysisTree root files              | 1. peaks_{source}.dat;</br> 2. peak_{source}.root, including uncalibrated spectrum with Gaussian fit on peaks;                                                                                                                                                                                        | 1. Only accept source type saved in "sources" folder</br> 2. Fit relative peaks (saved in "sources/") in uncalibrated spectrum. 2. It required "co60/co60_linfit.dat" exists; </br> 3. You don't have to run Step3, peaks_{source}.dat includes uncalibrated peak centroid and its related energy. |
| Step3          | Calibration.cxx           | 1. sources name                                                                        | 1. cal_pars.dat, including three arrays of calibration parameters; </br> 2. calibration.root, including TGraph with quad fit of each crystal and calibrated/uncalibrated summary plots.                                                                                                               | 1. It must run after "Calibration_HistMaker"</br> 2. It requires "peaks_{source}.dat" and "peaks_{source}.root" in "peaks/" folder;                                                                                                                                                                |


## sources
This folder includes sources peak information. Add "#" at the beginning of line that the energy you don't want to include in the calibration.

## macros
This folder so far only includes one python code to print out fwhm-relative info in co60\_linfit.dat. </br>
`python3 macros/checkfwhm.py co60_linfit.dat`
