import numpy as np

def read_resolution_file(filename):
  """
  Read resolution data file and extract 7 columns.

  Columns:
    0: Channel
    1: FWHM(Pu)
    2: FWHM(Am)
    3: FWHM(Cm)
    4: Resolution (%) at 5.8 MeV
    5: Gain
    6: Offset

  Returns:
    channel (int ndarray)
    fwhm_pu (ndarray)
    fwhm_am (ndarray)
    fwhm_cm (ndarray)
    res     (ndarray)
    gain    (ndarray)
    offset  (ndarray)
  """
  data = np.loadtxt(filename, comments="#")

  channel = data[:, 0].astype(int)
  fwhm_pu = data[:, 1]
  fwhm_am = data[:, 2]
  fwhm_cm = data[:, 3]
  res     = data[:, 4]
  gain    = data[:, 5]
  offset  = data[:, 6]

  return channel, fwhm_pu, fwhm_am, fwhm_cm, res, gain, offset

# ==================================================================== #
import matplotlib.pyplot as plt
import numpy as np

def compare_resolution(files):
  """
  Compare Res% channel-by-channel for multiple files.
  Channel == 0 will be removed before plotting.
  """

  plt.figure(figsize=(10, 6))

  for filename in files:
    channel, _, _, _, res, _, _ = read_resolution_file(filename)

    # --- remove channel == 0 ---
    mask = channel != 0
    channel = channel[mask]
    res     = res[mask]

    plt.plot(channel, res, marker='o', linestyle='-', label=filename)

  # --- x-axis zoom ---
  plt.xlim(960, 1071)

  # --- detector boundaries ---
  boundaries = {
    960:  "Detector 1 sector",
    991:  "Detector 1 ring",
    1015: "Detector 2 sector",
    1047: "Detector 2 ring",
  }

  for x, label in boundaries.items():
    plt.axvline(x=x, color='k', linestyle='--', linewidth=1)
    plt.text(
      x + 1,
      plt.ylim()[1] * 0.98,
      label,
      rotation=90,
      verticalalignment='top',
      fontsize=9
    )

  # --- grid: horizontal only ---
  plt.grid(axis='y')

  plt.xlabel("Channel")
  plt.ylabel("Resolution (%) at 5.8 MeV")
  plt.title("Channel-by-Channel Resolution Comparison")
  plt.legend()

  plt.tight_layout()
  plt.show()

# ==================================================================== #
if __name__ == "__main__":

  files = [
    "Res_Check_62337.dat",
    "Res_Check_62338.dat",
    "Res_Check_62341.dat",
  ]

  compare_resolution(files)

