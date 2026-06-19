#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_fluka_scored_dlet_composition_all8 as flcomp


TAG = os.environ.get("TAG", "FRAGLET_COMSCW_30M")
NCHUNKS = int(os.environ.get("NCHUNKS", "10"))

OUT = Path(f"plots_fragment_fluence_dose_profiles_{TAG}")
OUT.mkdir(exist_ok=True)

# Fragment-only components.
# Fluence forts:
#   41 DFLU_ZN, 42 TFLU_ZN, 43 H3FL_ZN, 46 H4FL_ZN, 48 LI6_ZN, 49 LI7_ZN
#
# Dose forts:
#   62 DODS_ZN, 63 TDOS_ZN, 64 H3DS_ZN, 65 H4DS_ZN, 66 LI6_DZ, 67 LI7_DZ
FRAGMENTS = [
    ("deuteron", 41, 62, "tab:blue"),
    ("triton",   42, 63, "tab:cyan"),
    ("helium3",  43, 64, "tab:green"),
    ("helium4",  46, 65, "tab:olive"),
    ("lithium6", 48, 66, "tab:orange"),
    ("lithium7", 49, 67, "tab:red"),
]


def read_sum_over_chunks(case, fort):
    cdir = Path(f"fluka_{case}_letmom")
    acc = np.zeros(flcomp.NZ)

    for chunk in range(1, NCHUNKS + 1):
        p = cdir / f"{case}_{TAG}_chunk{chunk}_fort{fort}.lis"
        acc += flcomp.read_usrbin_lis_values(p)

    # Each chunk is an independent FLUKA run normalized per primary.
    # Average the chunk profiles to recover the per-primary estimate.
    return acc / float(NCHUNKS)


def positive_ylim(curves):
    vals = []
    for y in curves:
        yy = np.asarray(y)
        yy = yy[np.isfinite(yy) & (yy > 0.0)]
        if yy.size:
            vals.append(yy)

    if not vals:
        return 1.0

    joined = np.concatenate(vals)
    return np.nanmax(joined) * 1.15


for case in flcomp.CASES:
    z = flcomp.z_axis_for_case(case)

    flu_curves = {}
    dose_curves = {}

    for name, flu_fort, dose_fort, color in FRAGMENTS:
        flu_curves[name] = read_sum_over_chunks(case, flu_fort)

        # FLUKA DOSE scorer is GeV/g/primary; convert to MeV/g/primary.
        dose_curves[name] = read_sum_over_chunks(case, dose_fort) * 1000.0

    # ------------------------------------------------------------------
    # 1) Fragment fluence
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 5))

    for name, flu_fort, dose_fort, color in FRAGMENTS:
        ax.plot(z, flu_curves[name], color=color, linewidth=1.6, label=name)

    ax.set_title(f"{case}\nFLUKA fragment fluence profiles, {TAG}")
    ax.set_xlabel("z [cm]")
    ax.set_ylabel("Fluence-like score [cm$^{-2}$ primary$^{-1}$]")
    ax.set_xlim(float(z[0]), float(z[-1]))
    ax.set_ylim(0.0, positive_ylim(flu_curves.values()))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()

    out_flu = OUT / f"{case}_fragment_fluence_profiles_{TAG}.png"
    fig.savefig(out_flu, dpi=200)
    plt.close(fig)

    # ------------------------------------------------------------------
    # 2) Fragment dose
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 5))

    for name, flu_fort, dose_fort, color in FRAGMENTS:
        ax.plot(z, dose_curves[name], color=color, linewidth=1.6, label=name)

    ax.set_title(f"{case}\nFLUKA fragment dose profiles, {TAG}")
    ax.set_xlabel("z [cm]")
    ax.set_ylabel("Dose [MeV/g/primary]")
    ax.set_xlim(float(z[0]), float(z[-1]))
    ax.set_ylim(0.0, positive_ylim(dose_curves.values()))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()

    out_dose = OUT / f"{case}_fragment_dose_profiles_{TAG}.png"
    fig.savefig(out_dose, dpi=200)
    plt.close(fig)

    print(case, "->", out_flu, "and", out_dose)

print("wrote", OUT)
