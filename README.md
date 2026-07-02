# macros_root

## Contents

1. [Overview](#overview)
   1. [Files](#files)
   2. [Fitting Function](#fitting-function)
   3. [Photopeak Area](#photopeak-area)
   4. [TSpectrum Background Removal](#tspectrum-background-removal)
2. [photopeakfit()](#photopeakfit)
   1. [Interface And Use](#interface-and-use)
   2. [Parameters, Initial Values, And Limits](#parameters-initial-values-and-limits)
   3. [Fit Mode And Options](#fit-mode-and-options)
   4. [PhotoPeakFit_Config.txt](#photopeakfit_configtxt)
3. [multipeakfit()](#multipeakfit)
   1. [Interface And Use](#interface-and-use-1)
   2. [Parameters, Initial Values, And Limits](#parameters-initial-values-and-limits-1)
   3. [Fit Mode And Options](#fit-mode-and-options-1)
   4. [MultiPeakFit_Config.txt](#multipeakfit_configtxt)

# Overview

`PhotoPeakFit.C` is a CERN ROOT macro for RadWare/GF3-style photopeak fitting. It can be loaded directly in ROOT without compilation and provides two public fitting calls:

- `photopeakfit()` for one photopeak.
- `multipeakfit()` for multiple photopeaks in one common fit range.

## Files

| File | Purpose |
| --- | --- |
| `PhotoPeakFit.C` | ROOT macro containing fit functions, config parsing, result printing, and drawing. |
| `Config/PhotoPeakFit_Config.txt` | Example config file for single-peak `photopeakfit()`. |
| `Config/MultiPeakFit_Config.txt` | Example config file for multi-peak `multipeakfit()`. |

## Fitting Function

The single-peak GF3-style function is

$$F(x) = G(x) + S(x) + Q(x) + B_{\mathrm{step}}(x),$$

where `G` is the Gaussian photopeak, `S` is the skew-tail photopeak, `Q` is the quadratic background, and `B_step` is the smoothed step background.

The definitions are

$$\sigma = \frac{W}{2.35482},$$

$$w = \frac{x - P}{\sqrt{2}\sigma},$$

$$y = \frac{W}{3.33021838\,\mathrm{BETA}},$$

$$G(x) = H\left(1 - \frac{R}{100}\right)e^{-w^2},$$

$$S(x) = H\frac{R}{100}\frac{\exp\left(\frac{x-P}{\mathrm{BETA}}\right)\mathrm{erfc}(w+y)}{\mathrm{erfc}(y)},$$

$$Q(x) = A + B\,x_c + C\,x_c^2,$$

$$B_{\mathrm{step}}(x) = H\,\mathrm{STEP}\,\frac{\mathrm{erfc}(w)}{200},$$

with

$$x_c = x - \frac{x_{\mathrm{low}} + x_{\mathrm{high}}}{2}.$$

For `multipeakfit()`, `A/B/C/R/BETA/STEP` are shared and each peak contributes its own `G_i`, `S_i`, and `B_step,i`:

$$F(x) = Q(x) + \sum_i \left[G_i(x) + S_i(x) + B_{\mathrm{step},i}(x)\right].$$

## Photopeak Area

The photopeak area is calculated from only the Gaussian and skew-tail photopeak terms. The quadratic background and smoothed step background are not included.

Define

$$d = \frac{\exp(-y^2)}{\mathrm{erfc}(y)}.$$

The continuous x-integrated area is

$$\mathrm{Area}_{x} = H\left[\frac{R}{100}\,2\,\mathrm{BETA}\,d + \left(1 - \frac{R}{100}\right)W\,1.06446705\right].$$

The printed area is converted to histogram-bin counts:

$$\mathrm{Area} = \frac{\mathrm{Area}_{x}}{\Delta x},$$

where `Delta x` is the bin width at the fitted centroid. The uncertainty is propagated from the ROOT covariance matrix. In `multipeakfit()`, each peak gets its own area and uncertainty from that peak component, not from slicing the summed curve.

## TSpectrum Background Removal

TSpectrum background removal is controlled by the config keyword `option`, independently from fit mode.

| Option | Behavior |
| --- | --- |
| `none` | Do not call `TSpectrum::Background`; fit the original histogram. |
| `global` | Estimate TSpectrum background over the full histogram x range, subtract it, then fit. |
| `local` | Estimate TSpectrum background over a local range, subtract it in the fitting range, then fit. |

If `option` is omitted, the default is `none`.

For `global`, the default iteration count is 20:

```txt
option = global
iteration = 20
```

For `local`, the default iteration count is:

```txt
iteration = clamp(round(2.5 * FWHM_bins), 6, 20)
```

The local range can be set explicitly:

```txt
option = local
range = 180,200
```

If `range` is omitted, the macro uses the fit range plus this padding:

```txt
max(5 * FWHM, 0.5 * fitWidth, iteration * binWidth)
```

# photopeakfit()

`photopeakfit()` fits one photopeak, prints the selected model and fit results, and draws the total fit plus background on the histogram.

Each call removes fit curves and legends previously drawn by this macro on the current canvas before drawing the new result.

## Interface And Use

Start ROOT in this directory, load the macro, open a ROOT file, draw the histogram, and fit one peak:

```cpp
.L PhotoPeakFit.C
TFile *f = TFile::Open("158_12_05_25_back_subtracted.root");
auto h = (TH1D*)f->Get("ParticleGates/Er/GammaEfficiency_Er");
h->Draw();
photopeakfit(h, 175, 210, 191.75);
```

Explicit mode:

```cpp
photopeakfit(h, 175, 210, 191.75, kPhotoPeakHighStat);
photopeakfit(h, 175, 210, 191.75, kPhotoPeakLowStat);
```

Config file:

```cpp
photopeakfit(h, 175, 210, 191.75, "Config/PhotoPeakFit_Config.txt");
photopeakfit(h, 175, 210, "Config/PhotoPeakFit_Config.txt");
```

Public interfaces:

```cpp
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 PhotoPeakFitMode mode = kPhotoPeakAuto);

int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 int mode);

int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 const char *configFile);

int photopeakfit(TH1 *hist, double fitLow, double fitHigh,
                 const char *configFile);
```

A C++ `PhotoPeakFitConfig` overload also exists for advanced macro control.

Inputs:

- `hist`: histogram pointer.
- `fitLow`, `fitHigh`: fit limits as x-values.
- `peak0`: initial peak position as an x-value.
- `mode`: optional fitting mode.
- `configFile`: optional text config file.

Position rule:

- If `peak0` is supplied in the function call, it defines the initial peak position.
- If `peak0` is omitted, `configFile` must define `init P`.
- `limit P` and `fix P` remain valid constraints in both cases.

## Parameters, Initial Values, And Limits

| Parameter | Meaning | Default initial value | Default limit |
| --- | --- | --- | --- |
| `A` | Constant background term. | `(y(fitLow) + y(fitHigh)) / 2` | none |
| `B` | Linear background slope in centered x. | `(y(fitHigh) - y(fitLow)) / (fitHigh - fitLow)` | none |
| `C` | Quadratic background coefficient. | `0` | none |
| `R` | Skew-tail fraction in percent. | `10` | `[0, 100]` |
| `BETA` | Skew-tail decay constant. | `W0 / 2` | `[1e-6, 10 * fitWidth]` |
| `STEP` | Smoothed step relative height. | `0.25` | `[0, 100]` |
| `P` | Peak centroid position. | `peak0`, or `init P` in config-only calls | `[fitLow, fitHigh]` |
| `W` | FWHM. | `sqrt(9 + 0.004 * peak0)` | `[1e-6, fitWidth]` |
| `H` | Peak height. | peak bin content minus estimated linear background | `[0, 10 * yMax]` |

## Fit Mode And Options

Fit modes:

| Mode | Behavior |
| --- | --- |
| `kPhotoPeakAuto` / `auto` | Try allowed candidates and choose the best reduced chi-square. |
| `kPhotoPeakHighStat` / `highstat` | Use the full 9-parameter function. |
| `kPhotoPeakLowStat` / `lowstat` | Use Gaussian plus linear background. |

Candidate model parts:

| Part | Parameters |
| --- | --- |
| Gaussian plus linear background | `A`, `B`, `P`, `W`, `H` |
| Skew tail | `R`, `BETA` |
| Smoothed step | `STEP` |
| Quadratic background | `C` |

`rootFitOption` is passed to ROOT:

```cpp
hist->Fit(func, rootFitOption);
```

The default `RQSN` means:

| ROOT option | Meaning |
| --- | --- |
| `R` | Use the TF1 fit range. |
| `Q` | Quiet fit output. |
| `S` | Return `TFitResultPtr` for covariance and status. |
| `N` | Do not store or draw ROOT's automatic fit function. |

`option` controls TSpectrum background removal and accepts `none`, `global`, or `local`.

## PhotoPeakFit_Config.txt

`Config/PhotoPeakFit_Config.txt` is the single-peak example config. It supports:

```txt
mode = auto
mode = highstat
mode = lowstat

rootFitOption = RQSN
option = none # global, local
iteration = 20
range = 180,200

init <par> = <value>
limit <par> = <low> <high>
fix <par> = <value>
```

Example:

```txt
mode = highstat
rootFitOption = RQSN
option = local
iteration = 12
range = 1320,1345
init P = 1332.5
limit P = 1330 1335
limit W = 0.1 20
fix R = 10
fix BETA = 1.5
```

Parameter names must be `A`, `B`, `C`, `R`, `BETA`, `STEP`, `P`, `W`, and `H`.

`init P` is required only for the config-only interface:

```cpp
photopeakfit(h, low, high, "Config/PhotoPeakFit_Config.txt");
```

When `peak0` is supplied in the function call, `peak0` wins and `init P` does not override it.

`mode` decides the active model. `init`, `limit`, and `fix` only apply to parameters active in that model. In `auto`, each candidate uses the applicable config lines; inactive config lines for the final model are reported as warnings.

# multipeakfit()

`multipeakfit()` fits multiple photopeaks in one histogram range. It shares `A/B/C/R/BETA/STEP` across all peaks, while each peak has its own `P[i]`, `W[i]`, and `H[i]` unless relative constraints or fixed parameters make `P[i]` or `W[i]` derived.

Each call removes fit curves, single-peak curves, background curves, TSpectrum curves, and legends previously drawn by this macro on the current canvas before drawing the new result.

## Interface And Use

Quick use:

```cpp
.L PhotoPeakFit.C
multipeakfit(h, 185, 205, {191.75, 194.20, 198.10});
multipeakfit(h, 185, 205, {191.75, 194.20, 198.10},
             "Config/MultiPeakFit_Config.txt");
multipeakfit(h, 185, 205, "Config/MultiPeakFit_Config.txt");
```

Public interfaces:

```cpp
int multipeakfit(TH1 *hist, double fitLow, double fitHigh,
                 std::initializer_list<double> peaks);

int multipeakfit(TH1 *hist, double fitLow, double fitHigh,
                 std::initializer_list<double> peaks,
                 const char *configFile);

int multipeakfit(TH1 *hist, double fitLow, double fitHigh,
                 const char *configFile);
```

A C++ `PhotoPeakFitConfig` overload and a C-style peak-array overload also exist for advanced macro control.

Inputs:

- `hist`: histogram pointer.
- `fitLow`, `fitHigh`: common fit limits as x-values.
- `peaks`: initial peak positions as x-values, in user-given order.
- `configFile`: optional multi-peak text config file.

Peak order is preserved. It is recommended to enter peaks from low x to high x. If the order is not increasing, the macro prints a warning but does not sort, because order defines indexed config entries such as `P[0]`, `W[0]`, and `H[0]`.

Position rule:

- If `peaks` are supplied in the function call, they define the initial/reference peak positions.
- If `peaks` are omitted, `configFile` must define contiguous `init P[0]`, `init P[1]`, ... values.
- `limit P[i]` and `fix P[i]` remain valid constraints in both cases.

## Parameters, Initial Values, And Limits

Global shared parameters:

| Parameter | Meaning | Default initial value | Default limit |
| --- | --- | --- | --- |
| `A` | Shared constant background term. | `(y(fitLow) + y(fitHigh)) / 2` | none |
| `B` | Shared linear background slope. | `(y(fitHigh) - y(fitLow)) / fitWidth` | none |
| `C` | Shared quadratic background coefficient. | `0` | none |
| `R` | Shared skew-tail fraction in percent. | `10` | `[0, 100]` |
| `BETA` | Shared skew-tail decay constant. | `W0[0] / 2` | `[1e-6, 10 * fitWidth]` |
| `STEP` | Shared step relative height. | `0.25` | `[0, 100]` |

Per-peak and relative parameters:

| Parameter | Meaning | Default initial value | Default limit |
| --- | --- | --- | --- |
| `P[i]` | Peak `i` centroid or reference centroid. | input peak `i`, or `init P[i]` in config-only calls | `[fitLow, fitHigh]` when fitted |
| `W[i]` | Peak `i` FWHM or reference FWHM. | `init W[i]`, else `sqrt(9 + 0.004 * P0_i)` | `[1e-6, fitWidth]` when fitted |
| `H[i]` | Peak `i` height. | peak bin content minus estimated linear background | `[0, 10 * yMax]` |
| `WSCALE` | Common FWHM scale when `relativeFwhm = true`. | `1.0` | `[1e-6, 10]`, or `limit WSCALE` |

When `relativePosition = false`, each non-fixed `P[i]` is fitted independently.

When `relativePosition = true`, positions are derived from one master position:

```txt
P_i = PMASTER + (P0_i - P0_0)
```

If any `fix P[i]` is set, the first fixed `P[i]` anchors all positions and no position parameter is fitted. This is the relative-position or position-scale behavior: peak spacings are fixed by function-input peaks, or by `init P[i]` in config-only calls.

When `relativeFwhm = false`, each non-fixed `W[i]` is fitted independently.

When `relativeFwhm = true`, non-fixed widths are derived from `WSCALE`:

```txt
W_i = WSCALE * W0_i
```

`fix W[i]` fixes only peak `i`; other non-fixed peaks still use `WSCALE`.

## Fit Mode And Options

The `mode`, candidate model parts, `rootFitOption`, and TSpectrum `option` keywords have the same meaning as in `photopeakfit()`.

Multi-peak defaults:

| Setting | Default |
| --- | --- |
| `mode` | `auto` |
| `rootFitOption` | `RQSN` |
| `option` | `none` |
| `relativePosition` | `false` |
| `relativeFwhm` | `true` |

`auto` chooses the shape/background complexity under the current relative constraints. It does not automatically choose whether positions or widths should be relative; those are user-controlled assumptions.

If `fix R = 0`, the skew tail is disabled and `BETA` is removed as an active fit parameter.

The printed output first lists fit metadata, then groups each peak's position, height, FWHM, and area with uncertainties. Raw TF1 parameters are printed after the grouped peak values.

## MultiPeakFit_Config.txt

`Config/MultiPeakFit_Config.txt` is a complete 3-peak example. After `relativePosition` and `relativeFwhm`, it is organized into three sections matching the single-peak config style:

- Initial values
- Limits
- Fixed parameters

`init P[i]` values are required only for the config-only interface:

```cpp
multipeakfit(h, low, high, "Config/MultiPeakFit_Config.txt");
```

When the peak list is supplied in the function call, the function-input positions win and `init P[i]` does not override them. The config-only interface requires contiguous `init P[0]`, `init P[1]`, ... values.

Supported syntax:

```txt
mode = auto
rootFitOption = RQSN
option = none # global, local
iteration = 20
range = 185,205

relativePosition = false
relativeFwhm = true

init WSCALE = 1.0
limit WSCALE = 0.5 2.0
fix WSCALE = 1.0

init P[0] = 191.75
limit P[0] = 190 193
fix P[0] = 191.75

init W[0] = 3.0
limit W[0] = 0.1 20
fix W[0] = 3.0

init H[0] = 800
limit H[0] = 0 1e9
fix H[0] = 800
```

In `relativePosition = true`, `limit P[0]` applies to `PMASTER` when there is no fixed position anchor. `limit P[i>0]` only applies when `relativePosition = false`.

In `relativeFwhm = true`, use `limit WSCALE` for the fitted width scale. `limit W[i]` only applies when `relativeFwhm = false`; `init W[i]` remains useful as the reference width `W0_i`.
