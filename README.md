# macros_root

## Contents

1. [PhotoPeakFit.C](#photopeakfitc)

## PhotoPeakFit.C

`PhotoPeakFit.C` is a CERN ROOT macro for fitting one RadWare/GF3-style photopeak in a histogram. It can be loaded directly in ROOT without compilation.

### Quick Use

Start ROOT in this directory, load the macro, open the ROOT file, draw the histogram, and fit one peak:

```cpp
.L PhotoPeakFit.C
TFile *f = TFile::Open("158_12_05_25_back_subtracted.root");
auto h = (TH1D*)f->Get("ParticleGates/Er/GammaEfficiency_Er");
h->Draw();
photopeakfit(h, 175, 210, 191.75);
```

The function signature is:

```cpp
int photopeakfit(TH1 *hist, double fitLow, double fitHigh, double peak0,
                 int mode = kPhotoPeakAuto)
```

Inputs:

- `hist`: pointer to the histogram to fit.
- `fitLow`: lower fitting limit as an x-value.
- `fitHigh`: upper fitting limit as an x-value.
- `peak0`: initial peak position as an x-value.
- `mode`: optional fit mode. The default is `kPhotoPeakAuto`.

The three numerical inputs are x-values, not bin numbers.

The macro prints the requested mode, selected fitting function, fit status, $\chi^2$, number of degrees of freedom, reduced $\chi^2$ when defined, photopeak area with uncertainty, and all fitted parameters with uncertainties. It draws the total fit as a red solid line and the background as a black dashed line on the original histogram plot.

### Fit Modes

The original 9-parameter fitting function is grouped into four parts:

1. `gaussian_linearBg`: normal Gaussian photopeak plus linear background. This is the required base function and has 5 free parameters: `Centroid`, `FWHM`, `Height`, `bg0`, and `bg1`.
2. `tail`: low-energy skew tail, controlled by `R` and `BETA`.
3. `step`: smoothed step function, controlled by `STEP`.
4. `quadBg`: quadratic background term, controlled by `bg2`.

Available public modes are:

| Mode | Behavior |
| --- | --- |
| `kPhotoPeakAuto` | Default. Starts from `gaussian_linearBg` and tries all free combinations of `tail`, `step`, and `quadBg`. |
| `kPhotoPeakHighStat` | Uses the full function: `gaussian_linearBg_tail_step_quadBg`. This is the original 9-parameter fit. |
| `kPhotoPeakLowStat` | Uses only `gaussian_linearBg`. |

The selected fitting function is printed after each fit. Example names are:

| Printed function | Included parts |
| --- | --- |
| `gaussian_linearBg` | Gaussian photopeak plus linear background |
| `gaussian_linearBg_tail` | Base function plus low-energy skew tail |
| `gaussian_linearBg_step_quadBg` | Base function plus step and quadratic background |
| `gaussian_linearBg_tail_step_quadBg` | Full 9-parameter function |

In `kPhotoPeakAuto`, candidates with optional parts are only tried when the fit range has more bins than free parameters, so the reduced $\chi^2$ is defined. The base `gaussian_linearBg` function is always allowed, even with only 4 or 5 bins, because it is the minimum useful photopeak model. Among candidates with a defined reduced $\chi^2$, auto mode chooses the fit with the smallest reduced $\chi^2$, preferring successful ROOT fit status when there is a successful candidate. If no candidate has a defined reduced $\chi^2$, auto mode falls back to `gaussian_linearBg`.

### Fitting Function

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

$$Q(x) = \mathrm{bg0} + \mathrm{bg1}\,x_c + \mathrm{bg2}\,x_c^2,$$

$$B_{\mathrm{step}}(x) = H\,\mathrm{STEP}\,\frac{\mathrm{erfc}(w)}{200}.$$

Here

$$x_c = x - \frac{x_{\mathrm{low}} + x_{\mathrm{high}}}{2}.$$

The background curve drawn by the macro is

$$B_{\mathrm{drawn}}(x) = Q(x) + B_{\mathrm{step}}(x).$$

### Parameters And Defaults

The parameter order is:

| Parameter | Meaning | Default initial value |
| --- | --- | --- |
| $\mathrm{bg0}$ | Constant part of the quadratic background | $\frac{y(x_{\mathrm{low}}) + y(x_{\mathrm{high}})}{2}$ |
| $\mathrm{bg1}$ | Linear background slope in centered x | $\frac{y(x_{\mathrm{high}}) - y(x_{\mathrm{low}})}{x_{\mathrm{high}} - x_{\mathrm{low}}}$ |
| $\mathrm{bg2}$ | Quadratic background coefficient | $0$ |
| $R$ | Skew-tail fraction in percent | $10$ |
| $\mathrm{BETA}$ | Skew-tail decay constant | $\frac{W_0}{2}$ |
| $\mathrm{STEP}$ | Smoothed step relative height | $0.25$ |
| $\mathrm{Centroid}$ | Peak centroid | `peak0` |
| $\mathrm{FWHM}$ | Full width at half maximum | $W_0 = \sqrt{9 + 0.004\,\mathrm{peak0}}$ |
| $\mathrm{Height}$ | Peak height | $y(\mathrm{peak0})$ minus estimated linear background |

Default parameter limits are:

| Parameter | Default limit |
| --- | --- |
| $\mathrm{bg0},\mathrm{bg1},\mathrm{bg2}$ | No explicit limit |
| $R$ | $[0,100]$ |
| $\mathrm{BETA}$ | $[10^{-6}, 10(x_{\mathrm{high}} - x_{\mathrm{low}})]$ |
| $\mathrm{STEP}$ | $[0,100]$ |
| $\mathrm{Centroid}$ | $[x_{\mathrm{low}}, x_{\mathrm{high}}]$ |
| $\mathrm{FWHM}$ | $[10^{-6}, x_{\mathrm{high}} - x_{\mathrm{low}}]$ |
| $\mathrm{Height}$ | $[0, 10\,y_{\max}]$ in the fitting range |

### Photopeak Area

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

where $V_{ij}$ is the covariance matrix returned by ROOT and $p_i$ are the fit parameters. Only derivatives with respect to $R$, $\mathrm{BETA}$, $W$, and $H$ are non-zero for this area formula.
