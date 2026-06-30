# macros_root

## Contents

1. [Overview](#overview)
   1. [Fitting Function](#fitting-function)
   2. [Parameters](#parameters)
   3. [Defaults And Limits](#defaults-and-limits)
   4. [Photopeak Area](#photopeak-area)
2. [photopeakfit()](#photopeakfit)
   1. [Files](#files)
   2. [Quick Use](#quick-use)
   3. [Function Interfaces](#function-interfaces)
   4. [Fit Modes](#fit-modes)
   5. [PhotoPeakFit_Config.txt](#photopeakfit_configtxt)

# Overview

`PhotoPeakFit.C` is a CERN ROOT macro for fitting one RadWare/GF3-style photopeak in a histogram. It can be loaded directly in ROOT without compilation.

## Fitting Function

The total fit function is

$$F(x) = G(x) + S(x) + Q(x) + B_{\mathrm{step}}(x),$$

where the four terms are the Gaussian photopeak, skew-tail photopeak, quadratic background, and smoothed step background.

The definitions are

$$\sigma = \frac{W}{2.35482},$$

$$w = \frac{x - P}{\sqrt{2}\sigma},$$

$$y = \frac{W}{3.33021838\,\mathrm{BETA}},$$

and

$$G(x) = H\left(1 - \frac{R}{100}\right)e^{-w^2},$$

$$S(x) = H\frac{R}{100}\frac{\exp\left(\frac{x-P}{\mathrm{BETA}}\right)\mathrm{erfc}(w+y)}{\mathrm{erfc}(y)},$$

$$Q(x) = A + B\,x_c + C\,x_c^2,$$

$$B_{\mathrm{step}}(x) = H\,\mathrm{STEP}\,\frac{\mathrm{erfc}(w)}{200}.$$

Here

$$x_c = x - \frac{x_{\mathrm{low}} + x_{\mathrm{high}}}{2}.$$

The background curve drawn by the macro is

$$B_{\mathrm{drawn}}(x) = Q(x) + B_{\mathrm{step}}(x).$$

## Parameters

Parameter names match the enum order in `PhotoPeakFit.C`:

| Parameter | Meaning |
| --- | --- |
| `A` | Constant term of the quadratic background. |
| `B` | Linear background slope in centered x. |
| `C` | Quadratic background coefficient. |
| `R` | Skew-tail fraction in percent. |
| `BETA` | Skew-tail decay constant. |
| `STEP` | Smoothed step relative height. |
| `P` | Peak centroid position. |
| `W` | FWHM, full width at half maximum. |
| `H` | Fitted peak height. |

## Defaults And Limits

Default initial values are:

| Parameter | Default initial value |
| --- | --- |
| `A` | $\frac{y(x_{\mathrm{low}}) + y(x_{\mathrm{high}})}{2}$ |
| `B` | $\frac{y(x_{\mathrm{high}}) - y(x_{\mathrm{low}})}{x_{\mathrm{high}} - x_{\mathrm{low}}}$ |
| `C` | $0$ |
| `R` | $10$ |
| `BETA` | $\frac{W_0}{2}$ |
| `STEP` | $0.25$ |
| `P` | `peak0` |
| `W` | $W_0 = \sqrt{9 + 0.004\,\mathrm{peak0}}$ |
| `H` | $y(\mathrm{peak0})$ minus estimated linear background |

Default parameter limits are:

| Parameter | Default limit |
| --- | --- |
| `A`, `B`, `C` | No explicit limit |
| `R` | $[0,100]$ |
| `BETA` | $[10^{-6}, 10(x_{\mathrm{high}} - x_{\mathrm{low}})]$ |
| `STEP` | $[0,100]$ |
| `P` | $[x_{\mathrm{low}}, x_{\mathrm{high}}]$ |
| `W` | $[10^{-6}, x_{\mathrm{high}} - x_{\mathrm{low}}]$ |
| `H` | $[0, 10\,y_{\max}]$ in the fitting range |

## Photopeak Area

The photopeak area is calculated from only the first two photopeak terms, $G(x)$ and $S(x)$. The quadratic background and smoothed step background are not included.

Define

$$d = \frac{\exp(-y^2)}{\mathrm{erfc}(y)}.$$

The continuous x-integrated area is

$$\mathrm{Area}_{x} = H\left[\frac{R}{100}\,2\,\mathrm{BETA}\,d + \left(1 - \frac{R}{100}\right)W\,1.06446705\right].$$

For the value printed by `photopeakfit`, this continuous area is converted to the histogram-bin-count convention:

$$\mathrm{Area} = \frac{\mathrm{Area}_{x}}{\Delta x},$$

where $\Delta x$ is the histogram bin width at the fitted centroid. This matches area estimates based on summing histogram bin contents. For a histogram with $\Delta x = 0.5$, the printed area is twice the continuous x-integrated area.

The uncertainty is propagated from the fit covariance matrix:

$$\sigma^2_{\mathrm{Area}} = \sum_i\sum_j \frac{\partial \mathrm{Area}}{\partial p_i} V_{ij} \frac{\partial \mathrm{Area}}{\partial p_j},$$

where $V_{ij}$ is the covariance matrix returned by ROOT and $p_i$ are the fit parameters. Only derivatives with respect to `R`, `BETA`, `W`, and `H` are non-zero for this area formula.

# photopeakfit()

`photopeakfit()` fits one photopeak, prints the selected model and fit results, and draws the total fit plus background on the histogram.

## Files

| File | Purpose |
| --- | --- |
| `PhotoPeakFit.C` | ROOT macro containing the fit function, mode selection, config parser, fit output, and drawing. |
| `Config/PhotoPeakFit_Config.txt` | Example user-editable config file for mode, initial values, limits, and fixed parameters. |

## Quick Use

Start ROOT in this directory, load the macro, open the ROOT file, draw the histogram, and fit one peak:

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
```

When a config file is used, `photopeakfit` reads the file each time it is called. You can keep the ROOT session open, edit `Config/PhotoPeakFit_Config.txt`, and run the same command again.

The three numerical inputs are x-values, not bin numbers.

## Function Interfaces

The public interfaces are:

```cpp
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 PhotoPeakFitMode mode = kPhotoPeakAuto);

int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 const char *configFile);

int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 const PhotoPeakFitConfig &config);
```

Inputs:

- `hist`: pointer to the histogram to fit.
- `fitLow`: lower fitting limit as an x-value.
- `fitHigh`: upper fitting limit as an x-value.
- `peak0`: initial peak position as an x-value.
- `mode`: optional fit mode. The default is `kPhotoPeakAuto`.
- `configFile`: optional text file containing mode, initial values, limits, and fixed parameters.
- `config`: optional C++ configuration object for macro-level control.

The macro prints the requested mode, selected fitting function, fit status, $\chi^2$, number of degrees of freedom, reduced $\chi^2$ when defined, photopeak area with uncertainty, and all fitted parameters with uncertainties. It draws the total fit as a red solid line and the background as a black dashed line on the original histogram plot.

## Fit Modes

The macro supports three public fit modes:

| Mode | Behavior |
| --- | --- |
| `kPhotoPeakAuto` | Default. Tries allowed candidate models and selects the best reduced $\chi^2$. |
| `kPhotoPeakHighStat` | Uses the full 9-parameter function. |
| `kPhotoPeakLowStat` | Uses only the Gaussian plus linear background base function. |

The 9-parameter fitting function is grouped into four parts:

1. `gaussian_linearBg`: normal Gaussian photopeak plus linear background. This is the required base function and has 5 free parameters: `P`, `W`, `H`, `A`, and `B`.
2. `tail`: low-energy skew tail, controlled by `R` and `BETA`.
3. `step`: smoothed step function, controlled by `STEP`.
4. `quadBg`: quadratic background term, controlled by `C`.

Available candidate names printed by the macro include:

| Printed function | Included parts |
| --- | --- |
| `gaussian_linearBg` | Gaussian photopeak plus linear background. |
| `gaussian_linearBg_tail` | Base function plus low-energy skew tail. |
| `gaussian_linearBg_step_quadBg` | Base function plus step and quadratic background. |
| `gaussian_linearBg_tail_step_quadBg` | Full 9-parameter function. |

In `kPhotoPeakAuto`, candidates with optional parts are only tried when the fit range has more bins than free parameters, so the reduced $\chi^2$ is defined. User fixed parameters reduce the free-parameter count for this check. The base `gaussian_linearBg` function is always allowed, even with only 4 or 5 bins, because it is the minimum useful photopeak model.

Among candidates with a defined reduced $\chi^2$, auto mode chooses the fit with the smallest reduced $\chi^2$, preferring successful ROOT fit status when there is a successful candidate. If no candidate has a defined reduced $\chi^2$, auto mode falls back to `gaussian_linearBg`.

## PhotoPeakFit_Config.txt

`Config/PhotoPeakFit_Config.txt` is an example config file. Blank lines and lines starting with `#` are ignored.

Minimal config:

```txt
mode = auto
```

Example with user constraints:

```txt
mode = highstat
init P = 1332.5
limit P = 1330 1335
limit W = 0.1 20
fix R = 10
fix BETA = 1.5
```

Supported syntax:

```txt
mode = auto
mode = highstat
mode = lowstat

init <par> = <value>
limit <par> = <low> <high>
fix <par> = <value>
```

Parameter names must match `PhotoPeakFit.C`: `A`, `B`, `C`, `R`, `BETA`, `STEP`, `P`, `W`, and `H`.

Conflict rule:

`mode` decides the active model. `init`, `limit`, and `fix` only apply to parameters active in that model. For example, if `mode = lowstat` and `fix R = 10`, `R` is inactive because the low-stat model has no skew tail. The fit remains low-stat, and the `fix R` line is ignored with a warning.

In `mode = auto`, each candidate model only uses config lines for parameters active in that candidate. After the final model is selected, inactive config lines for that model are reported as warnings.
