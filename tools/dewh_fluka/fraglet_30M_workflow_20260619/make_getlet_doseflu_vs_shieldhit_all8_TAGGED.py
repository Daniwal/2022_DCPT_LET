
#!/usr/bin/env python3
"""
Build FLUKA-vs-SHIELD-HIT comparison CSV files for all 8 GETLET benchmark cases.

Purpose
-------
This script is the central post-processing script for the current trusted
GETLET comparison workflow. It reads already-converted FLUKA USRBIN `.lis`
files and SHIELD-HIT depth-profile `.dat` files, puts them on the same depth
axis, reconstructs LET quantities from FLUKA LET-moment scorers, and writes
one comparison CSV per case plus one all-case summary CSV.

The resulting per-case CSVs are used by:

    plot_all_getlet_vs_shieldhit.py

to produce the trusted 48-plot package:

    GETLET_TRUSTED_48_PLOTS_20260603

Inputs
------
For each case, the script expects a FLUKA case directory:

    /home/dewh/runs/let_benchmark/fluka_<case>_letmom/

containing:

    <case>_letmom_step2_prod.inp
        FLUKA production input. Used only to recover the Z_narrow depth grid.

    getlet_doseflu_prod_fort<N>.lis
        ASCII-converted FLUKA USRBIN outputs, produced with FLUKA `usbrea`.

The script also expects SHIELD-HIT results in:

    /home/dewh/2022_DCPT_LET/data/sh12a/results/<case>/

containing the relevant NB_Z_narrow dose, fluence, and LET reference files.

Main physical quantities
------------------------
FLUKA quantities:
    EDEP_all
        All-particle energy deposition from FLUKA USRBIN ENERGY scoring.
        The raw FLUKA DOSE value is treated as GeV/g/primary.

    PDEP_protons
        Proton energy deposition from custom/proton-filtered scoring.
        Also treated as GeV/g/primary.

    AFLU_all
        All-particle fluence-like USRBIN score.

    PHI_protons
        Proton fluence moment denominator.

    PRI_primary
        Primary-proton fluence moment denominator, using LTRACK primary filtering.

    PHL1, PHL2, PWL1, PWL2
        Proton LET moment scorers. The suffixes represent first and second
        LET-weighted moments in local material or water-reference stopping power.

Reconstructed LET definitions
-----------------------------
The script reconstructs track-averaged and dose-averaged LET from moments:

    TLET = sum(phi_i * L_i) / sum(phi_i)

    DLET = sum(phi_i * L_i^2) / sum(phi_i * L_i)

where `phi_i` represents the appropriate fluence/track contribution and `L_i`
is the LET/stopping-power value used in the custom FLUKA scoring routine.

For example:

    FL_TLET_water = PWL1_water / PHI_protons
    FL_DLET_water = PWL2_water / PWL1_water

and similarly for primary protons using PRI/PWR scorers.

Units
-----
Depth:
    cm

Dose:
    FLUKA raw DOSE/PDOSE values are converted from:

        GeV/g/primary

    to:

        MeV/g/primary

    using the solid-water density:

        rho_SOLWATR = 1.032 g/cm^3

    Conversion:

        dose[MeV/g/primary] = EDEP[GeV/g/primary] * 1000 / rho

LET:
    keV/um

Fluence:
    cm^-2 primary^-1, where appropriate for USRBIN fluence scoring.

Masking and metrics
-------------------
The script defines a dose mask:

    FL_DOSE_all_norm > 0.05

i.e. bins above 5% of the maximum FLUKA all-particle dose.

This mask is used for LET and comparison metrics to avoid interpreting
low-dose distal-tail bins where LET ratios can become numerically unstable.

Important assumptions and caveats
---------------------------------
1. The FLUKA `.lis` files must contain the physical data matrix first.
   This script reads values after the first "Data follow" marker and assumes
   the first 205 numeric values correspond to the Z_narrow depth bins.

2. The script assumes 205 depth bins for the GETLET Z_narrow scoring region.

3. The FLUKA dose conversion uses solid-water density 1.032 g/cm^3.
   This is correct for the current SOLWATR/SWF benchmark comparison, but should
   be revisited if the geometry material changes.

4. LET quantities are only as correct as the custom FLUKA moment scorers that
   produced the corresponding fort files.

5. This script does not run FLUKA and does not run SHIELD-HIT. It only combines
   already-produced results.

Outputs
-------
For each case:

    fluka_<case>_letmom/<case>_getlet_doseflu_vs_SHIELDHIT.csv

Global summary:

    /home/dewh/runs/let_benchmark/getlet_doseflu_vs_SHIELDHIT_all8_summary.csv
"""
import re
import csv
import os
from pathlib import Path
import numpy as np
# Root directory for this FLUKA/SHIELD-HIT benchmark work on Grendel.
# All generated FLUKA case folders live below this directory, e.g.
#   /home/dewh/runs/let_benchmark/fluka_plan02_field01_geoD_mono_letmom
ROOT = Path("/home/dewh/runs/let_benchmark")

TAG = os.environ.get("TAG", "FRAGLET_COMSCW_30M")
INPUT_TAG = os.environ.get("INPUT_TAG", TAG)
NCHUNKS = int(os.environ.get("NCHUNKS", "10"))

# Root directory for the reference SHIELD-HIT12A results from the project repo.
# Each case has a subfolder here containing NB_Z_narrow_*.dat files.
SHROOT = Path("/home/dewh/2022_DCPT_LET/data/sh12a/results")

# Conversion factor for FLUKA DOSE scorer 228:
#
#   FLUKA DOSE [GeV/g/primary] -> [MeV/g/primary]
#
# Explanation:
#   scorer 228 = DOSE = energy deposited per unit mass, GeV/g.
#   Therefore only the energy-unit conversion is needed:
#
#       1 GeV/g = 1000 MeV/g
#
# No density division is applied here. The older workflow used scorer 208
# DOSE, which is treated as GeV/g and converted directly to MeV/g by multiplying by 1000
# density. That is no longer correct once fort.21/fort.27 use scorer 228.
GEV_PER_G_TO_MEV_PER_G = 1000.0

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

# Mapping between physical scorer names and the FLUKA fort numbers produced by
# the current trusted GETLET custom executable.
#
# The corresponding ASCII files are expected to be named:
#   getlet_doseflu_prod_fort<N>.lis
#
# The fort numbers come from the order of USRBIN scoring cards in the FLUKA
# input file. Do not change this mapping unless the scoring-card order or output
# units are intentionally changed and revalidated.
#
# Moment naming convention:
#   PHI / PHL / PWL:
#       all protons, i.e. primary + secondary protons
#
#   PRI / PRL / PWR:
#       primary protons only, selected in the custom FLUSCW routine using
#       LTRACK .EQ. 1
#
#   local:
#       LET evaluated in the material of the current region/step
#
#   water:
#       LET evaluated using water-reference stopping power
#
#   1:
#       first LET moment, sum(phi_i * L_i)
#
#   2:
#       second LET moment, sum(phi_i * L_i^2)
FORT = {
    "EDEP_all": 21,
    "PHI_protons": 22,
    "PHL1_local": 23,
    "PHL2_local": 24,
    "PWL1_water": 25,
    "PWL2_water": 26,
    "PDEP_protons": 27,
    "AFLU_all": 28,
    "PRI_primary": 29,
    "PRL1_primary_local": 30,
    "PRL2_primary_local": 31,
    "PWR1_primary_water": 32,
    "PWR2_primary_water": 33,
    "PRDO_primary": 34,
}

def values_after_data(path, n=205):
    """
    Extract the physical USRBIN data values from a FLUKA `usbrea` ASCII file.

    FLUKA binary scoring files are converted to human-readable `.lis` files
    using `usbrea`. A typical `.lis` file contains a text header followed by
    the actual scored data matrix after a line containing:

        Data follow

    For this project, each Z_narrow scorer has:

        nx = 1
        ny = 1
        nz = 205

    so the physical data matrix should contain 205 bin values.

    Parameters
    ----------
    path:
        Path to the `.lis` file.

    n:
        Number of physical depth bins expected. For the GETLET Z_narrow
        benchmark this is 205.

    Returns
    -------
    numpy.ndarray
        Array of length `n` containing the first physical USRBIN data matrix.

    Notes
    -----
    The function deliberately reads only the first `n` numbers after the
    physical-data marker. This avoids accidentally using later auxiliary
    matrices, such as percentage-error matrices, if present.
    """

    # Store all numeric values found after the "Data follow" marker.
    vals = []

    # This flag becomes True once we have reached the physical data matrix.
    take = False

    # Read the file as text. `errors="ignore"` avoids failures from odd
    # non-UTF characters that can appear in files generated by legacy tools.
    txt = Path(path).read_text(errors="ignore").splitlines()

    for line in txt:
        # The first "Data follow" line marks the start of the physical scored
        # USRBIN matrix.
        if "Data follow" in line:
            take = True
            continue

        # Ignore everything before the physical data matrix begins.
        if not take:
            continue

        # These are descriptive lines printed by `usbrea`, not numeric data.
        if "accurate deposition" in line or "this is" in line:
            continue

        # Extract Fortran/scientific-notation values such as:
        #   1.2345E-06
        #
        # The regex intentionally matches the numeric matrix values while
        # ignoring text labels and headers.
        vals.extend(float(x) for x in re.findall(r"[-+]?\d+\.\d+E[-+]\d+", line, re.I))

    # A valid Z_narrow profile must contain at least 205 values.
    if len(vals) < n:
        raise RuntimeError(f"Only found {len(vals)} values in {path}")

    # Return exactly the expected number of depth bins.
    return np.array(vals[:n], dtype=float)

def z_grid_from_input(inp):
    """
    Reconstruct the Z_narrow depth-bin centers from the FLUKA input file.

    The FLUKA input contains the DOSE_ZN USRBIN scoring card. In FLUKA USRBIN
    syntax, the first card gives the upper corner of the scoring mesh, and the
    continuation card gives the lower corner and number of bins.

    For this project, the relevant cards look schematically like:

        USRBIN ... x_max y_max z_max DOSE_ZN
        USRBIN ... x_min y_min z_min n_x n_y n_z &

    This function reads z_min, z_max, and n_z, then returns the bin centers:

        z_i = z_min + (i + 0.5) * dz

    where:

        dz = (z_max - z_min) / n_z

    Using the FLUKA input as the source of the depth grid avoids hard-coding
    different phantom lengths for the mono, SOBP, and ramp cases.
    """

    # Read the FLUKA production input as plain text.
    lines = Path(inp).read_text().splitlines()

    for i, line in enumerate(lines):
        # Find the USRBIN card defining the all-particle energy deposition
        # scorer in the Z_narrow central scoring region.
        if "DOSE_ZN" in line and line.startswith("USRBIN"):
            # `top` is the first USRBIN card: contains the upper z boundary.
            top = line.split()

            # `bot` is the continuation card immediately below: contains the
            # lower z boundary and number of z bins.
            bot = lines[i + 1].split()

            # FLUKA fixed/free-format USRBIN layout in this generated input:
            #   top[6] = z_max
            #   bot[3] = z_min
            #   bot[6] = n_z
            zmax = float(top[6])
            zmin = float(bot[3])
            nz = int(float(bot[6]))

            # Convert bin edges to bin-center coordinates.
            dz = (zmax - zmin) / nz
            return zmin + (np.arange(nz) + 0.5) * dz

    # If this fails, the input may not be a trusted GETLET Z_narrow input.
    raise RuntimeError(f"Could not find DOSE_ZN grid in {inp}")

def load_sh(path, z_target, scale=1.0):
    """
    Load a SHIELD-HIT depth-profile file and interpolate it onto the FLUKA grid.

    SHIELD-HIT `.dat` files used here are expected to have at least two columns:

        column 0: depth z [cm]
        column 1: scored quantity

    The FLUKA and SHIELD-HIT depth grids are intended to match, but interpolation
    is used to make the comparison robust against tiny floating-point formatting
    differences in the depth coordinates.

    Parameters
    ----------
    path:
        Path to the SHIELD-HIT `.dat` file.

    z_target:
        Depth-bin centers from the FLUKA Z_narrow grid.

    scale:
        Optional multiplicative scale factor applied to the SHIELD-HIT scored
        value before interpolation. The default is 1.0.

    Returns
    -------
    numpy.ndarray
        SHIELD-HIT quantity evaluated on the FLUKA depth grid.
    """

    arr = np.loadtxt(path)

    # Interpolate the SHIELD-HIT scored value, column 1, onto the FLUKA depth
    # grid. Values outside the SHIELD-HIT range should not occur for properly
    # matched cases.
    return np.interp(z_target, arr[:, 0], arr[:, 1] * scale)


def norm(v):
    """
    Normalize an array to its maximum finite value.

    This is used for shape comparisons, e.g. normalized dose and normalized
    fluence. It should not be used when absolute scaling is physically relevant.

    If the maximum is zero or non-finite, return an array of NaNs instead of
    silently dividing by an invalid value.
    """

    m = np.nanmax(v)

    if m == 0 or not np.isfinite(m):
        return np.full_like(v, np.nan)

    return v / m


def safe_div(a, b):
    """
    Divide two arrays while avoiding division by zero.

    LET moment ratios can have zero denominators in low-fluence or low-dose
    bins, especially in the distal tail. Rather than producing infinities, this
    function returns NaN where the denominator is zero.

    Examples in this script:
        TLET = first LET moment / fluence moment
        DLET = second LET moment / first LET moment
    """

    # Start with NaN everywhere; fill only bins with non-zero denominator.
    out = np.full_like(a, np.nan, dtype=float)

    mask = b != 0
    out[mask] = a[mask] / b[mask]

    return out

def metrics(a, b, mask):
    """
    Compute simple difference metrics between two depth profiles.

    The comparison is evaluated only in bins selected by `mask`. In this script,
    that is usually the region where the FLUKA all-particle dose is above 5% of
    its maximum value. This avoids letting low-dose distal-tail bins dominate
    the comparison.

    Parameters
    ----------
    a:
        First profile, usually FLUKA.

    b:
        Second profile, usually SHIELD-HIT.

    mask:
        Boolean array selecting the depth bins included in the metric.

    Returns
    -------
    dict
        MAD:
            Mean absolute difference.

        RMSD:
            Root-mean-square difference.

        MAXABS:
            Maximum absolute difference.

    Notes
    -----
    The metrics are absolute differences in the units of the input profiles.
    For normalized dose/fluence profiles, the values are dimensionless.
    For absolute dose, the values are in MeV/g/primary.
    For LET, the values are in keV/um.
    """

    # Difference profile in the selected region.
    d = a[mask] - b[mask]

    # Remove NaN/inf values caused by empty bins or safe divisions.
    d = d[np.isfinite(d)]

    return {
        "MAD": float(np.mean(np.abs(d))) if len(d) else np.nan,
        "RMSD": float(np.sqrt(np.mean(d * d))) if len(d) else np.nan,
        "MAXABS": float(np.max(np.abs(d))) if len(d) else np.nan,
    }


# Accumulate one dictionary of summary metrics per benchmark case.
# These rows are written at the end to:
#
#   getlet_doseflu_vs_SHIELDHIT_all8_summary.csv
summary_rows = []


for case in CASES:
    # FLUKA case directory for this benchmark case.
    # Example:
    #   /home/dewh/runs/let_benchmark/fluka_plan02_field01_geoD_mono_letmom
    cdir = ROOT / f"fluka_{case}_letmom"

    # SHIELD-HIT reference result directory for the same case.
    # Example:
    #   /home/dewh/2022_DCPT_LET/data/sh12a/results/plan02_field01_geoD_mono
    shdir = SHROOT / case

    # Recover the FLUKA Z_narrow depth grid directly from the production input.
    # This avoids hard-coding the depth range, which differs between mono/SOBP
    # and ramp-style cases.
    z = z_grid_from_input(cdir / f"{case}_letmom_step2_prod.inp")

    # Dictionary holding all FLUKA depth profiles for this case.
    # Keys are physical scorer names from the FORT mapping above.
    fl = {}

    for key, fort in FORT.items():
        # Read and sum chunked ASCII-converted FLUKA USRBIN files from the
        # tagged multi-chunk production run, e.g.
        #
        #   <case>_<TAG>_chunk<chunk>_fort<N>.lis
        #
        # For TAG=FRAGLET_COMSCW_30M and NCHUNKS=10 this gives the full
        # 30M-primary-per-case result.
        acc = None

        for chunk in range(1, NCHUNKS + 1):
            lis_path = cdir / f"{case}_{INPUT_TAG}_chunk{chunk}_fort{fort}.lis"
            vals = values_after_data(lis_path)

            if acc is None:
                acc = vals.copy()
            else:
                acc += vals

        # Average independent per-primary FLUKA chunks to recover the full-run per-primary estimate.
        fl[key] = acc / float(NCHUNKS)

    # Convert FLUKA DOSE USRBIN scorer 228 from GeV/g/primary to MeV/g/primary.
    #
    # fort.21 = all-particle dose
    # fort.27 = proton dose, selected with AUXSCORE PROTON
    fl["DOSE_all_MeV_g"] = fl["EDEP_all"] * GEV_PER_G_TO_MEV_PER_G
    fl["DOSE_protons_MeV_g"] = fl["PDEP_protons"] * GEV_PER_G_TO_MEV_PER_G
    fl["DOSE_primary_MeV_g"] = fl["PRDO_primary"] * GEV_PER_G_TO_MEV_PER_G

    # -------------------------------------------------------------------------
    # Reconstruct FLUKA LET quantities from custom LET-moment scorers.
    # -------------------------------------------------------------------------
    #
    # The custom FLUSCW scoring routine accumulates LET moments rather than
    # writing LET directly.
    #
    # Track-averaged LET:
    #
    #   TLET = sum(phi_i * L_i) / sum(phi_i)
    #
    # Dose-averaged LET:
    #
    #   DLET = sum(phi_i * L_i^2) / sum(phi_i * L_i)
    #
    # The `safe_div` helper returns NaN in bins where the denominator is zero,
    # which avoids artificial infinities in low-fluence distal-tail bins.

    # All-proton LET in local material.
    # Includes primary + secondary protons.
    fl["TLET_local"] = safe_div(fl["PHL1_local"], fl["PHI_protons"])
    fl["DLET_local"] = safe_div(fl["PHL2_local"], fl["PHL1_local"])

    # All-proton LET using water-reference stopping power.
    # These are the main quantities compared to SHIELD-HIT water LET.
    fl["TLET_water"] = safe_div(fl["PWL1_water"], fl["PHI_protons"])
    fl["DLET_water"] = safe_div(fl["PWL2_water"], fl["PWL1_water"])

    # Primary-proton LET in local material.
    # Primary protons are selected in the custom Fortran routine by LTRACK .EQ. 1.
    fl["TLET_primary_local"] = safe_div(fl["PRL1_primary_local"], fl["PRI_primary"])
    fl["DLET_primary_local"] = safe_div(fl["PRL2_primary_local"], fl["PRL1_primary_local"])

    # Primary-proton LET using water-reference stopping power.
    fl["TLET_primary_water"] = safe_div(fl["PWR1_primary_water"], fl["PRI_primary"])
    fl["DLET_primary_water"] = safe_div(fl["PWR2_primary_water"], fl["PWR1_primary_water"])

    # -------------------------------------------------------------------------
    # Load SHIELD-HIT reference dose and fluence profiles.
    # -------------------------------------------------------------------------
    #
    # These files come from the SHIELD-HIT Z_narrow output:
    #
    #   NB_Z_narrow_dose_p1.dat = fluence, all particles
    #   NB_Z_narrow_dose_p2.dat = dose, all particles
    #   NB_Z_narrow_dose_p3.dat = dose, protons
    #   NB_Z_narrow_dose_p4.dat = fluence, primary protons
    #   NB_Z_narrow_dose_p5.dat = fluence, all protons
    #
    # The SHIELD-HIT profiles are interpolated onto the FLUKA z grid so that
    # all arrays below have identical length and matching depth coordinates.
    sh = {}

    sh["fluence_all"] = load_sh(shdir / "NB_Z_narrow_dose_p1.dat", z)
    sh["dose_all_MeV_g"] = load_sh(shdir / "NB_Z_narrow_dose_p2.dat", z)
    sh["dose_protons_MeV_g"] = load_sh(shdir / "NB_Z_narrow_dose_p3.dat", z)
    sh["fluence_primary"] = load_sh(shdir / "NB_Z_narrow_dose_p4.dat", z)
    sh["fluence_protons"] = load_sh(shdir / "NB_Z_narrow_dose_p5.dat", z)

    # Normalized SHIELD-HIT profiles for shape-only comparisons.
    # These should not be interpreted as absolute dose or absolute fluence.
    sh["fluence_all_norm"] = norm(sh["fluence_all"])
    sh["dose_all_norm"] = norm(sh["dose_all_MeV_g"])
    sh["dose_protons_norm"] = norm(sh["dose_protons_MeV_g"])
    sh["fluence_primary_norm"] = norm(sh["fluence_primary"])
    sh["fluence_protons_norm"] = norm(sh["fluence_protons"])

    # -------------------------------------------------------------------------
    # Load SHIELD-HIT water-reference LET profiles.
    # -------------------------------------------------------------------------
    #
    # The relevant SHIELD-HIT files are:
    #
    #   NB_Z_narrow_LET_water_p6.dat = TLET, all protons, water reference
    #   NB_Z_narrow_LET_water_p3.dat = DLET, all protons, water reference
    #   NB_Z_narrow_LET_water_p5.dat = TLET, primary protons, water reference
    #   NB_Z_narrow_LET_water_p2.dat = DLET, primary protons, water reference
    #
    # The values are multiplied by 0.1 to convert the SHIELD-HIT output units
    # used here into keV/um, matching the FLUKA GETLET moment reconstruction.
    sh["TLET_water"] = load_sh(shdir / "NB_Z_narrow_LET_water_p6.dat", z, scale=0.1)
    sh["DLET_water"] = load_sh(shdir / "NB_Z_narrow_LET_water_p3.dat", z, scale=0.1)
    sh["TLET_primary_water"] = load_sh(shdir / "NB_Z_narrow_LET_water_p5.dat", z, scale=0.1)
    sh["DLET_primary_water"] = load_sh(shdir / "NB_Z_narrow_LET_water_p2.dat", z, scale=0.1)

    # -------------------------------------------------------------------------
    # Build normalized FLUKA profiles for shape-only comparisons.
    # -------------------------------------------------------------------------
    #
    # Absolute quantities are kept separately above. These normalized arrays are
    # used only for plotting/metrics where the profile shape is the relevant
    # comparison, not the absolute scale.
    fl_norm = {
        "EDEP_all_norm": norm(fl["EDEP_all"]),
        "PDEP_protons_norm": norm(fl["PDEP_protons"]),
        "PRDO_primary_norm": norm(fl["PRDO_primary"]),
        "AFLU_all_norm": norm(fl["AFLU_all"]),
        "PHI_protons_norm": norm(fl["PHI_protons"]),
        "PRI_primary_norm": norm(fl["PRI_primary"]),
    }

    # -------------------------------------------------------------------------
    # Define the analysis mask.
    # -------------------------------------------------------------------------
    #
    # LET values can become numerically unstable in very low-dose distal-tail
    # bins because the fluence or first-moment denominators become very small.
    #
    # We therefore use the same trusted convention as in the plotting workflow:
    # keep only bins where the normalized FLUKA all-particle dose is above 5%.
    #
    # This mask is stored in the per-case CSV so plotting scripts and future
    # audits can see exactly which bins were included in LET metrics.
    dose_mask = fl_norm["EDEP_all_norm"] > 0.05

    # -------------------------------------------------------------------------
    # Assemble the per-depth output table.
    # -------------------------------------------------------------------------
    #
    # `out_cols` and `header` must stay in exactly the same order.
    #
    # Each output row corresponds to one Z_narrow depth bin. The columns include:
    #
    #   - depth
    #   - FLUKA absolute and normalized dose
    #   - SHIELD-HIT absolute and normalized dose
    #   - FLUKA and SHIELD-HIT fluence diagnostics
    #   - FLUKA reconstructed LET quantities
    #   - SHIELD-HIT reference LET quantities
    #   - the 5% FLUKA-dose mask
    #
    # This CSV is intentionally broad: it contains enough information to make
    # dose, fluence, LET, and diagnostic plots without re-reading the raw `.lis`
    # or SHIELD-HIT `.dat` files.
    out_cols = [
        z,

        # All-particle dose.
        fl["EDEP_all"], fl["DOSE_all_MeV_g"], sh["dose_all_MeV_g"],
        fl_norm["EDEP_all_norm"], sh["dose_all_norm"],

        # Proton dose.
        fl["PDEP_protons"], fl["DOSE_protons_MeV_g"], sh["dose_protons_MeV_g"],
        fl_norm["PDEP_protons_norm"], sh["dose_protons_norm"],

        # Primary-proton dose.
        fl["PRDO_primary"], fl["DOSE_primary_MeV_g"], fl_norm["PRDO_primary_norm"],

        # All-particle fluence diagnostics.
        fl["AFLU_all"], sh["fluence_all"],
        fl_norm["AFLU_all_norm"], sh["fluence_all_norm"],

        # All-proton fluence diagnostics.
        fl["PHI_protons"], sh["fluence_protons"],
        fl_norm["PHI_protons_norm"], sh["fluence_protons_norm"],

        # Primary-proton fluence diagnostics.
        fl["PRI_primary"], sh["fluence_primary"],
        fl_norm["PRI_primary_norm"], sh["fluence_primary_norm"],

        # FLUKA local-material LET.
        fl["TLET_local"], fl["DLET_local"],

        # All-proton water-reference LET: FLUKA vs SHIELD-HIT.
        fl["TLET_water"], sh["TLET_water"],
        fl["DLET_water"], sh["DLET_water"],

        # Primary-proton LET.
        fl["TLET_primary_local"], fl["DLET_primary_local"],
        fl["TLET_primary_water"], sh["TLET_primary_water"],
        fl["DLET_primary_water"], sh["DLET_primary_water"],

        # Store mask as 0/1 integer for CSV readability.
        dose_mask.astype(int),
    ]

    # Column names for the per-case CSV.
    #
    # Keep this list synchronized with `out_cols` above. A mismatch would make
    # downstream plotting scripts read the wrong physical quantity.
    header = ",".join([
        "depth_cm",
        "FL_DOSE_all_GeV_g", "FL_dose_all_MeV_g", "SH_dose_all_MeV_g",
        "FL_EDEP_all_norm", "SH_dose_all_norm",
        "FL_PDOSE_protons_GeV_g", "FL_dose_protons_MeV_g", "SH_dose_protons_MeV_g",
        "FL_PDEP_protons_norm", "SH_dose_protons_norm",
        "FL_PRDO_primary_GeV_g", "FL_dose_primary_MeV_g", "FL_PRDO_primary_norm",
        "FL_AFLU_all", "SH_fluence_all", "FL_AFLU_all_norm", "SH_fluence_all_norm",
        "FL_PHI_protons", "SH_fluence_protons", "FL_PHI_protons_norm", "SH_fluence_protons_norm",
        "FL_PRI_primary", "SH_fluence_primary", "FL_PRI_primary_norm", "SH_fluence_primary_norm",
        "FL_TLET_local", "FL_DLET_local",
        "FL_TLET_water", "SH_TLET_water",
        "FL_DLET_water", "SH_DLET_water",
        "FL_TLET_primary_local", "FL_DLET_primary_local",
        "FL_TLET_primary_water", "SH_TLET_primary_water",
        "FL_DLET_primary_water", "SH_DLET_primary_water",
        "mask_5pct_FL_DOSE",
    ])

    # Combine the one-dimensional depth profiles into a 2D table:
    #
    #   rows    = depth bins
    #   columns = quantities listed in `header`
    out = np.column_stack(out_cols)

    # Write the per-case comparison CSV into the corresponding FLUKA case folder.
    out_csv = cdir / f"{case}_getlet_doseflu_vs_SHIELDHIT_{TAG}.csv"
    np.savetxt(out_csv, out, delimiter=",", header=header, comments="")

    # -------------------------------------------------------------------------
    # Build one summary-metric row for this case.
    # -------------------------------------------------------------------------
    #
    # The detailed per-depth CSV above is the main data product. This summary row
    # is a compact case-level overview used to quickly check agreement across
    # all 8 benchmark cases.
    row = {"case": case}

    # Dictionary of FLUKA-vs-SHIELD-HIT comparisons.
    #
    # Each entry has:
    #
    #   name: (FLUKA_profile, SHIELD-HIT_profile, mask)
    #
    # The `metrics` helper then computes:
    #
    #   MAD, RMSD, MAXABS
    #
    # inside the masked dose region.
    comparisons = {
        # Absolute dose comparisons in MeV/g/primary.
        "dose_all_abs_MeV_g": (fl["DOSE_all_MeV_g"], sh["dose_all_MeV_g"], dose_mask),
        "dose_protons_abs_MeV_g": (fl["DOSE_protons_MeV_g"], sh["dose_protons_MeV_g"], dose_mask),

        # Shape-only dose comparisons after normalization to each profile maximum.
        "dose_all_shape": (fl_norm["EDEP_all_norm"], sh["dose_all_norm"], dose_mask),
        "dose_protons_shape": (fl_norm["PDEP_protons_norm"], sh["dose_protons_norm"], dose_mask),

        # Shape-only fluence diagnostics.
        "fluence_all_shape": (fl_norm["AFLU_all_norm"], sh["fluence_all_norm"], dose_mask),
        "fluence_protons_shape": (fl_norm["PHI_protons_norm"], sh["fluence_protons_norm"], dose_mask),
        "fluence_primary_shape": (fl_norm["PRI_primary_norm"], sh["fluence_primary_norm"], dose_mask),

        # LET comparisons in keV/um.
        "TLET_water": (fl["TLET_water"], sh["TLET_water"], dose_mask),
        "DLET_water": (fl["DLET_water"], sh["DLET_water"], dose_mask),
        "TLET_primary_water": (fl["TLET_primary_water"], sh["TLET_primary_water"], dose_mask),
        "DLET_primary_water": (fl["DLET_primary_water"], sh["DLET_primary_water"], dose_mask),
    }

    # Compute and store all summary metrics.
    for name, (a, b, m) in comparisons.items():
        mm = metrics(a, b, m)
        row[f"{name}_MAD"] = mm["MAD"]
        row[f"{name}_RMSD"] = mm["RMSD"]
        row[f"{name}_MAXABS"] = mm["MAXABS"]

    # Store maximum LET values inside the same 5% dose mask.
    #
    # This gives a compact check of LET magnitude without being dominated by
    # numerical spikes in low-dose bins outside the clinically/physically useful
    # comparison region.
    row["FL_TLET_water_max"] = float(np.nanmax(fl["TLET_water"][dose_mask]))
    row["FL_DLET_water_max"] = float(np.nanmax(fl["DLET_water"][dose_mask]))
    row["SH_TLET_water_max"] = float(np.nanmax(sh["TLET_water"][dose_mask]))
    row["SH_DLET_water_max"] = float(np.nanmax(sh["DLET_water"][dose_mask]))

    row["FL_TLET_primary_water_max"] = float(np.nanmax(fl["TLET_primary_water"][dose_mask]))
    row["FL_DLET_primary_water_max"] = float(np.nanmax(fl["DLET_primary_water"][dose_mask]))
    row["SH_TLET_primary_water_max"] = float(np.nanmax(sh["TLET_primary_water"][dose_mask]))
    row["SH_DLET_primary_water_max"] = float(np.nanmax(sh["DLET_primary_water"][dose_mask]))

    # Add this case to the global all-case summary table.
    summary_rows.append(row)

    # Print a short progress line with the most important current diagnostics.
    print(
        f"{case}: wrote {out_csv.name} "
        f"TLETw_MAD={row['TLET_water_MAD']:.6g} "
        f"DLETw_MAD={row['DLET_water_MAD']:.6g} "
        f"dose_abs_MAD={row['dose_all_abs_MeV_g_MAD']:.6g} "
        f"dose_shape_MAD={row['dose_all_shape_MAD']:.6g}"
    )


# -----------------------------------------------------------------------------
# Write the all-case summary CSV.
# -----------------------------------------------------------------------------
#
# This file gives one row per benchmark case and is useful for quickly checking
# whether any case is anomalous after changing the custom FLUKA routine,
# conversion scripts, or plotting workflow.
summary_csv = ROOT / f"getlet_doseflu_vs_SHIELDHIT_all8_summary_{TAG}.csv"

# Use the keys from the first row so column order follows the construction above.
fieldnames = list(summary_rows[0].keys())

with open(summary_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(summary_rows)

print()
print("Wrote", summary_csv)