#!/usr/bin/env python3
import sys
import statistics as stats

# ANSI colors for red + bold
RED_BOLD = "\033[1;31m"
RESET = "\033[0m"

if len(sys.argv) != 2:
    print("Usage: python3 check_fwhm.py co60_linfit.dat")
    sys.exit(1)

filename = sys.argv[1]

good_values = []

print("=== FWHM Report ===")

with open(filename, "r") as f:
    for line in f:
        line = line.strip()

        # Skip empty or commented lines
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < 2:
            continue

        array_num = parts[0]

        # Extract the second last column (FWHM)
        try:
            fwhm = float(parts[-2])
        except ValueError:
            continue

        # Check if this value is bad (> 5 keV)
        if fwhm > 5.0:
            # Print in red + bold
            print(f"{RED_BOLD}ArrayNumber:\t{array_num}\tFWHM(1332):\t{fwhm:.6f} keV{RESET}")
        else:
            good_values.append(fwhm)
            print(f"ArrayNumber:\t{array_num}\tFWHM(1332):\t{fwhm:.6f} keV")

# If no good values found
if not good_values:
    print("\nNo good FWHM values found (<= 5 keV).")
    sys.exit(0)

# Compute statistics
mean_val = stats.mean(good_values)
std_val = stats.stdev(good_values) if len(good_values) > 1 else 0.0
min_val = min(good_values)
max_val = max(good_values)

# Final summary output
print("\n=== RESULTS (good FWHM values only) ===")
print(f"Cores = {len(good_values)}")
print(f"Mean  = {mean_val:.6f} keV")
print(f"Std   = {std_val:.6f}")
print(f"Min   = {min_val:.6f} keV")
print(f"Max   = {max_val:.6f} keV")

