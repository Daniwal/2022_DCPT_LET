#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_fluka_scored_dlet_composition_all8 as flcomp


flcomp.TAG = os.environ.get("TAG", flcomp.TAG)
flcomp.NCHUNKS = int(os.environ.get("NCHUNKS", flcomp.NCHUNKS))

TAG = flcomp.TAG
NCHUNKS = flcomp.NCHUNKS

SHROOT = Path(os.environ.get("SHROOT", "/home/dewh/2022_DCPT_LET/data/sh12a/results"))
OUT = Path("plots_fluka_scored_dlet_composition_vs_shieldhit_LOCAL_{}".format(TAG))
COMPONENT_COLORS = {
    "proton":   "lightgray",
    "deuteron": "tab:blue",
    "triton":   "tab:cyan",
    "helium3":  "tab:green",
    "helium4":  "tab:olive",
    "lithium6": "tab:orange",
    "lithium7": "tab:red",
}
OUT.mkdir(exist_ok=True)


def load_sh(path, scale=0.1):
    a = np.loadtxt(path)
    return a[:, 0], a[:, 1] * scale


def interpolate_to_fluka(z_sh, y_sh, z_fl):
    return np.interp(z_fl, z_sh, y_sh, left=np.nan, right=np.nan)


def ylimit_from(*arrays):
    vals = []
    for a in arrays:
        finite = np.asarray(a)[np.isfinite(a)]
        finite = finite[finite >= 0.0]
        if finite.size:
            vals.append(finite)

    if not vals:
        return 5.0

    joined = np.concatenate(vals)
    return max(5.0, np.nanpercentile(joined, 99.5) * 1.15)


for case in flcomp.CASES:
    cdir = Path("fluka_{}_letmom".format(case))
    z = flcomp.z_axis_for_case(case)

    let1_by_component = {}
    let2_by_component = {}

    for name, fort_let1, fort_let2 in flcomp.COMPONENTS:
        let1_sum = np.zeros(flcomp.NZ)
        let2_sum = np.zeros(flcomp.NZ)

        for chunk in range(1, NCHUNKS + 1):
            f1 = cdir / "{}_{}_chunk{}_fort{}.lis".format(case, TAG, chunk, fort_let1)
            f2 = cdir / "{}_{}_chunk{}_fort{}.lis".format(case, TAG, chunk, fort_let2)

            let1_sum += flcomp.read_usrbin_lis_values(f1)
            let2_sum += flcomp.read_usrbin_lis_values(f2)

        let1_by_component[name] = let1_sum
        let2_by_component[name] = let2_sum

    denom = np.zeros(flcomp.NZ)
    total_let2 = np.zeros(flcomp.NZ)

    for name, _, _ in flcomp.COMPONENTS:
        denom += let1_by_component[name]
        total_let2 += let2_by_component[name]

    fluka_scored_dlet = flcomp.safe_divide(total_let2, denom)

    shdir = SHROOT / case

    z_sh, sh_dlet_all_local = load_sh(shdir / "NB_Z_narrow_LET_p1.dat", scale=0.1)
    _, sh_dlet_proton_local = load_sh(shdir / "NB_Z_narrow_LET_p3.dat", scale=0.1)
    z_dose, sh_dose = load_sh(shdir / "NB_Z_narrow_dose_p2.dat", scale=1.0)

    sh_dlet_all_i = interpolate_to_fluka(z_sh, sh_dlet_all_local, z)
    sh_dlet_proton_i = interpolate_to_fluka(z_sh, sh_dlet_proton_local, z)

    sh_dose_i = interpolate_to_fluka(z_dose, sh_dose, z)
    dose_mask = sh_dose_i / np.nanmax(sh_dose_i) > 0.05
    denom_mask = denom > 0.0
    mask = dose_mask & denom_mask

    component_curves = []
    component_labels = []
    component_colors = []

    for name, _, _ in flcomp.COMPONENTS:
        comp = flcomp.safe_divide(let2_by_component[name], denom)
        comp = np.where(mask, comp, 0.0)
        component_curves.append(comp)
        component_labels.append(name)
        component_colors.append(COMPONENT_COLORS[name])

    fluka_line = np.where(mask, fluka_scored_dlet, np.nan)
    sh_all_line = np.where(dose_mask, sh_dlet_all_i, np.nan)
    sh_proton_line = np.where(dose_mask, sh_dlet_proton_i, np.nan)

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.stackplot(
        z,
        component_curves,
        labels=component_labels,
        colors=component_colors,
        alpha=0.85,
    )
    ax.plot(z, fluka_line, linewidth=1.6, label="FLUKA scored sum")
    ax.plot(z, sh_all_line, "--", linewidth=1.8, label="SHIELD-HIT all charged, local")
    ax.plot(z, sh_proton_line, ":", linewidth=1.8, label="SHIELD-HIT protons, local")

    ax.set_title("{}\nDLET composition: FLUKA scored species vs SHIELD-HIT local".format(case))
    ax.set_xlabel("z [cm]")
    ax.set_ylabel("DLET contribution [keV/um]")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax.set_xlim(float(z[0]), float(z[-1]))
    ax.set_ylim(0.0, ylimit_from(fluka_line, sh_all_line, sh_proton_line))

    fig.tight_layout()

    out = OUT / "{}_FLUKA_scored_DLET_composition_vs_SHIELDHIT_local.png".format(case)
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print("{} -> {}".format(case, out))

print("wrote", OUT)
