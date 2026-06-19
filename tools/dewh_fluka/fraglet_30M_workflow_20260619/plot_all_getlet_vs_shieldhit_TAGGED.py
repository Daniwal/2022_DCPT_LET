#!/usr/bin/env python3
"""
Plot FLUKA-vs-SHIELD-HIT dose, fluence, and LET comparisons for all 8 GETLET cases.

Purpose
-------
This script is the plotting companion to:

    make_getlet_doseflu_vs_shieldhit_all8.py

It reads the per-case comparison CSVs produced by that script and generates
one set of diagnostic/comparison PNGs per benchmark case.

Current trusted output set
--------------------------
The current trusted plotting workflow produces:

    8 cases × 6 plots per case = 48 PNG files

in:

    /home/dewh/runs/let_benchmark/plots_getlet_vs_shieldhit_20260528

The output folder name still contains the older date `20260528`, but the
trusted 48-plot package was updated and archived on 2026-06-03.

Inputs
------
For each case, the script expects:

    /home/dewh/runs/let_benchmark/fluka_<case>_letmom/<case>_getlet_doseflu_vs_SHIELDHIT.csv

These CSVs contain already-aligned FLUKA and SHIELD-HIT depth profiles,
including dose, fluence diagnostics, LET, and a 5% FLUKA-dose mask.

Outputs
-------
For each case, the script writes:

    01_dose_absolute.png
        Absolute FLUKA vs SHIELD-HIT dose comparison.

    02_fluence.png
        Normalized fluence diagnostic comparison.

    03_fluence_absolute.png
        Absolute fluence diagnostic comparison.

    04_LET_water_protons_masked.png
        All-proton TLET/DLET in water reference, masked to >5% FLUKA dose.

    05_LET_local_vs_water_FLUKA_masked.png
        FLUKA local-material LET compared against FLUKA water-reference LET,
        masked to >5% FLUKA dose.

    06_LET_water_primary_masked.png
        Primary-proton TLET/DLET in water reference, masked to >5% FLUKA dose.

Plotting conventions
--------------------
FLUKA:
    solid lines

SHIELD-HIT:
    dashed lines

Matching physical quantities:
    same color across FLUKA and SHIELD-HIT

LET plots:
    only show bins where mask_5pct_FL_DOSE is true

Important caveats
-----------------
1. This script does not recompute any physics quantities. It only plots columns
   already written to the comparison CSVs.

2. The LET masking is intentional. Unmasked low-dose distal-tail bins can
   produce visually misleading LET spikes from tiny denominators.

3. The normalized fluence plot is diagnostic and shape-based. The absolute
   fluence plot should be used when absolute scale matters.

4. If the CSV column names change, this plotting script must be updated.
"""
import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Root directory for the benchmark workflow on Grendel.
ROOT = Path("/home/dewh/runs/let_benchmark")

TAG = os.environ.get("TAG", "FRAGLET_COMSCW_30M")

# Output folder for the GETLET comparison plots.
#
# Historical note:
# The folder name contains 20260528 because it was created earlier, but the
# current trusted 48-plot set was regenerated and archived on 2026-06-03.
OUT = ROOT / f"plots_getlet_vs_shieldhit_{TAG}"
OUT.mkdir(exist_ok=True)

# Plot colors.
#
# Convention:
#   - same physical quantity gets the same color across FLUKA and SHIELD-HIT
#   - code identity is shown by linestyle:
#       FLUKA      = solid
#       SHIELD-HIT = dashed
COL_DOSE_ALL = "tab:blue"
COL_DOSE_PROTONS = "tab:orange"

COL_FLUENCE_ALL = "tab:blue"
COL_FLUENCE_PROTONS = "tab:orange"
COL_FLUENCE_PRIMARY = "tab:green"

COL_TLET = "tab:blue"
COL_DLET = "tab:orange"

CASES = [
    "plan01_field01_geoA_SOBPcent",
    "plan01_field01_geoB_SOBP95",
    "plan01_field01_geoC_SOBP74",
    "plan02_field01_geoD_mono",
    "plan03_field01_geoA_rampFull",
    "plan03_field02_geoA_rampFull",
    "plan04_field01_geoA_rampMiddle",
    "plan04_field02_geoA_rampMiddle",
]

def load_case(case):
    """
    Load one per-case FLUKA-vs-SHIELD-HIT comparison CSV.

    The CSV is produced by:

        make_getlet_doseflu_vs_shieldhit_all8.py

    It contains one row per depth bin and one column per scored/reconstructed
    quantity. The first line is a comma-separated header.
    """

    p = ROOT / f"fluka_{case}_letmom" / f"{case}_getlet_doseflu_vs_SHIELDHIT_{TAG}.csv"

    # Read column names from the first line.
    header = p.read_text().splitlines()[0].split(",")

    # Read numeric data from all remaining lines.
    data = np.loadtxt(p, delimiter=",", skiprows=1)

    return header, data


def c(header, name):
    """
    Return the column index for a named quantity in the comparison CSV.

    Keeping column lookup name-based is safer than hard-coding integer column
    positions, especially while the CSV format is still evolving.
    """

    return header.index(name)

for case in CASES:
    # Load the per-depth comparison table for this benchmark case.
    #
    # h = list of CSV column names
    # a = numeric data array with shape:
    #
    #     rows    = Z_narrow depth bins
    #     columns = quantities listed in h
    h, a = load_case(case)

    # Depth-bin centers in cm.
    z = a[:, c(h, "depth_cm")]

    # Boolean mask selecting the trusted LET plotting region:
    #
    #   FLUKA all-particle dose > 5% of its maximum
    #
    # LET plots use this mask to avoid low-dose distal-tail bins where moment
    # ratios can become visually misleading.
    mask = a[:, c(h, "mask_5pct_FL_DOSE")].astype(bool)

    # -------------------------------------------------------------------------
    # 1) Absolute dose comparison
    # -------------------------------------------------------------------------
    plt.figure(figsize=(8,5))
    plt.plot(z, a[:, c(h, "FL_dose_all_MeV_g")],
             color=COL_DOSE_ALL, linestyle="-",
             label="FLUKA dose all")
    plt.plot(z, a[:, c(h, "SH_dose_all_MeV_g")],
             color=COL_DOSE_ALL, linestyle="--",
             label="SHIELD-HIT dose all")
    plt.plot(z, a[:, c(h, "FL_dose_protons_MeV_g")],
             color=COL_DOSE_PROTONS, linestyle="-",
             label="FLUKA dose protons")
    plt.plot(z, a[:, c(h, "SH_dose_protons_MeV_g")],
             color=COL_DOSE_PROTONS, linestyle="--",
             label="SHIELD-HIT dose protons")
    plt.xlabel("Depth z [cm]")
    plt.ylabel("Dose [MeV/g/primary]")
    plt.title(f"{case}: absolute dose comparison")
    plt.grid(True, which="major", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT / f"{case}_01_dose_absolute.png", dpi=180)
    plt.close()

    # -------------------------------------------------------------------------
    # 2) Normalized fluence diagnostic comparison
    # -------------------------------------------------------------------------
    #
    # This plot compares the shapes of the fluence-like quantities after each
    # curve has been normalized to its own maximum. It is useful for checking
    # spatial agreement, but not absolute scaling.
    plt.figure(figsize=(8,5))
    plt.plot(z, a[:, c(h, "FL_AFLU_all_norm")],
             color=COL_FLUENCE_ALL, linestyle="-",
             label="FLUKA fluence all: AFLU_ZN")
    plt.plot(z, a[:, c(h, "SH_fluence_all_norm")],
             color=COL_FLUENCE_ALL, linestyle="--",
             label="SHIELD-HIT fluence all")
    plt.plot(z, a[:, c(h, "FL_PHI_protons_norm")],
             color=COL_FLUENCE_PROTONS, linestyle="-",
             label="FLUKA fluence protons: PHI_ZN")
    plt.plot(z, a[:, c(h, "SH_fluence_protons_norm")],
             color=COL_FLUENCE_PROTONS, linestyle="--",
             label="SHIELD-HIT fluence protons")
    plt.plot(z, a[:, c(h, "FL_PRI_primary_norm")],
             color=COL_FLUENCE_PRIMARY, linestyle="-",
             label="FLUKA fluence primary: PRI_ZN")
    plt.plot(z, a[:, c(h, "SH_fluence_primary_norm")],
             color=COL_FLUENCE_PRIMARY, linestyle="--",
             label="SHIELD-HIT fluence primary")
    plt.xlabel("Depth z [cm]")
    plt.ylabel("Normalized fluence-like quantity")
    plt.title(f"{case}: fluence diagnostics")
    plt.grid(True, which="major", alpha=0.25)    
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT / f"{case}_02_fluence.png", dpi=180)
    plt.close()

    # -------------------------------------------------------------------------
    # 3) Absolute fluence diagnostic comparison
    # -------------------------------------------------------------------------
    #
    # This plot keeps the absolute fluence scale. It is the appropriate plot to
    # inspect when checking whether FLUKA and SHIELD-HIT agree in magnitude.
    plt.figure(figsize=(8,5))
    plt.plot(z, a[:, c(h, "FL_AFLU_all")],
             color=COL_FLUENCE_ALL, linestyle="-",
             label="FLUKA fluence all")
    plt.plot(z, a[:, c(h, "SH_fluence_all")],
             color=COL_FLUENCE_ALL, linestyle="--",
             label="SHIELD-HIT fluence all")
    plt.plot(z, a[:, c(h, "FL_PHI_protons")],
             color=COL_FLUENCE_PROTONS, linestyle="-",
             label="FLUKA fluence protons")
    plt.plot(z, a[:, c(h, "SH_fluence_protons")],
             color=COL_FLUENCE_PROTONS, linestyle="--",
             label="SHIELD-HIT fluence protons")
    plt.plot(z, a[:, c(h, "FL_PRI_primary")],
             color=COL_FLUENCE_PRIMARY, linestyle="-",
             label="FLUKA fluence primary")
    plt.plot(z, a[:, c(h, "SH_fluence_primary")],
             color=COL_FLUENCE_PRIMARY, linestyle="--",
             label="SHIELD-HIT fluence primary")
    plt.xlabel("Depth z [cm]")
    plt.ylabel("Fluence [cm$^{-2}$ primary$^{-1}$]")
    plt.title(f"{case}: absolute fluence comparison")
    plt.grid(True, which="major", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT / f"{case}_03_fluence_absolute.png", dpi=180)
    plt.close()

     # -------------------------------------------------------------------------
    # 4) All-proton LET in water reference, masked to >5% FLUKA dose
    # -------------------------------------------------------------------------
    #
    # This is the main LET comparison plot for the current GETLET workflow.
    #
    # It compares:
    #   FLUKA reconstructed TLET/DLET for all protons
    #   SHIELD-HIT TLET/DLET for all protons
    #
    # "All protons" means primary + secondary protons, but not deuterons,
    # tritons, helium fragments, or heavier ions.
    #
    # The plotted points are masked to the trusted >5% FLUKA-dose region.
    plt.figure(figsize=(8,5))
    plt.plot(z[mask], a[mask, c(h, "FL_TLET_water")],
             color=COL_TLET, linestyle="-",
             label="FLUKA TLET water")
    plt.plot(z[mask], a[mask, c(h, "SH_TLET_water")],
             color=COL_TLET, linestyle="--",
             label="SHIELD-HIT TLET water")
    plt.plot(z[mask], a[mask, c(h, "FL_DLET_water")],
             color=COL_DLET, linestyle="-",
             label="FLUKA DLET water")
    plt.plot(z[mask], a[mask, c(h, "SH_DLET_water")],
             color=COL_DLET, linestyle="--",
             label="SHIELD-HIT DLET water")
    plt.xlabel("Depth z [cm]")
    plt.ylabel("LET [keV/µm]")
    plt.title(f"{case}: proton LET, water reference, >5% dose")
    plt.grid(True, which="major", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT / f"{case}_04_LET_water_protons_masked.png", dpi=180)
    plt.close()

    # -------------------------------------------------------------------------
    # 5) FLUKA local-material LET vs water-reference LET, masked
    # -------------------------------------------------------------------------
    #
    # This is a FLUKA-only diagnostic plot.
    #
    # It shows how much the reconstructed LET changes when the stopping power is
    # evaluated in the local transport material instead of using water-reference
    # stopping power. This matters because SHIELD-HIT LET_water is water
    # referenced, while the FLUKA geometry includes solid water and PMMA.
    plt.figure(figsize=(8,5))
    plt.plot(z[mask], a[mask, c(h, "FL_TLET_local")],
             color=COL_TLET, linestyle="-",
             label="FLUKA TLET local material")
    plt.plot(z[mask], a[mask, c(h, "FL_DLET_local")],
             color=COL_DLET, linestyle="-",
             label="FLUKA DLET local material")
    plt.plot(z[mask], a[mask, c(h, "FL_TLET_water")],
             color=COL_TLET, linestyle="--",
             label="FLUKA TLET water reference")
    plt.plot(z[mask], a[mask, c(h, "FL_DLET_water")],
             color=COL_DLET, linestyle="--",
             label="FLUKA DLET water reference")
    plt.xlabel("Depth z [cm]")
    plt.ylabel("LET [keV/µm]")
    plt.title(f"{case}: FLUKA local-material vs water-reference LET >5% dose")
    plt.legend(fontsize=8)
    plt.grid(True, which="major", alpha=0.25)
    plt.tight_layout()
    plt.savefig(OUT / f"{case}_05_LET_local_vs_water_FLUKA_masked.png", dpi=180)
    plt.close()

    # -------------------------------------------------------------------------
    # 6) Primary-proton LET in water reference, masked to >5% FLUKA dose
    # -------------------------------------------------------------------------
    #
    # This plot compares primary-proton-only LET between FLUKA and SHIELD-HIT.
    #
    # In the FLUKA custom routine, primary protons are selected using:
    #
    #   LTRACK .EQ. 1
    #
    # This plot is useful for separating primary-beam LET behavior from the
    # contribution of secondary protons.
    plt.figure(figsize=(8,5))
    plt.plot(z[mask], a[mask, c(h, "FL_TLET_primary_water")],
             color=COL_TLET, linestyle="-",
             label="FLUKA TLET primary water")
    plt.plot(z[mask], a[mask, c(h, "SH_TLET_primary_water")],
             color=COL_TLET, linestyle="--",
             label="SHIELD-HIT TLET primary water")
    plt.plot(z[mask], a[mask, c(h, "FL_DLET_primary_water")],
             color=COL_DLET, linestyle="-",
             label="FLUKA DLET primary water")
    plt.plot(z[mask], a[mask, c(h, "SH_DLET_primary_water")],
             color=COL_DLET, linestyle="--",
             label="SHIELD-HIT DLET primary water")
    plt.xlabel("Depth z [cm]")
    plt.ylabel("LET [keV/µm]")
    plt.title(f"{case}: primary-proton LET, water reference, >5% dose")
    plt.grid(True, which="major", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(OUT / f"{case}_06_LET_water_primary_masked.png", dpi=180)
    plt.close()

print(f"Wrote plots to: {OUT}")
