#!/usr/bin/env python3

import os
import csv
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_fluka_scored_dlet_composition_all8 as flcomp


TAG = os.environ.get("TAG", "FRAGLET_COMSCW_1M")
NCHUNKS = int(os.environ.get("NCHUNKS", "10"))

SHROOT = Path(os.environ.get("SHROOT", "/home/dewh/2022_DCPT_LET/data/sh12a/results"))
OUT = Path("plots_fluka_cumulative_scored_LET_vs_shieldhit_LOCAL_{}".format(TAG))
OUT.mkdir(exist_ok=True)

COMPONENTS = [
    ("proton",   22, 23, 24),
    ("deuteron", 41, 52, 56),
    ("triton",   42, 53, 57),
    ("helium3",  43, 54, 58),
    ("helium4",  46, 55, 59),
    ("lithium6", 48, 50, 60),
    ("lithium7", 49, 51, 61),
]


def load_sh(path, scale=0.1):
    a = np.loadtxt(path)
    return a[:, 0], a[:, 1] * scale


def interp_to_fluka(z_sh, y_sh, z_fl):
    return np.interp(z_fl, z_sh, y_sh, left=np.nan, right=np.nan)


def safe_divide(num, den):
    out = np.full_like(num, np.nan, dtype=float)
    mask = den > 0.0
    out[mask] = num[mask] / den[mask]
    return out


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


rows = []

for case in flcomp.CASES:
    cdir = Path("fluka_{}_letmom".format(case))
    z = flcomp.z_axis_for_case(case)

    flu_total = np.zeros(flcomp.NZ)
    let1_total = np.zeros(flcomp.NZ)
    let2_total = np.zeros(flcomp.NZ)

    component_flu_sums = {}
    component_let1_sums = {}
    component_let2_sums = {}

    for name, fort_flu, fort_let1, fort_let2 in COMPONENTS:
        flu = np.zeros(flcomp.NZ)
        let1 = np.zeros(flcomp.NZ)
        let2 = np.zeros(flcomp.NZ)

        for chunk in range(1, NCHUNKS + 1):
            f_flu = cdir / "{}_{}_chunk{}_fort{}.lis".format(case, TAG, chunk, fort_flu)
            f_l1 = cdir / "{}_{}_chunk{}_fort{}.lis".format(case, TAG, chunk, fort_let1)
            f_l2 = cdir / "{}_{}_chunk{}_fort{}.lis".format(case, TAG, chunk, fort_let2)

            flu += flcomp.read_usrbin_lis_values(f_flu)
            let1 += flcomp.read_usrbin_lis_values(f_l1)
            let2 += flcomp.read_usrbin_lis_values(f_l2)

        flu_total += flu
        let1_total += let1
        let2_total += let2

        component_flu_sums[name] = float(np.sum(flu))
        component_let1_sums[name] = float(np.sum(let1))
        component_let2_sums[name] = float(np.sum(let2))

    fluka_tlet = safe_divide(let1_total, flu_total)
    fluka_dlet = safe_divide(let2_total, let1_total)

    shdir = SHROOT / case

    z_sh, sh_dlet_all_local = load_sh(shdir / "NB_Z_narrow_LET_p1.dat", scale=0.1)
    _, sh_dlet_proton_local = load_sh(shdir / "NB_Z_narrow_LET_p3.dat", scale=0.1)
    _, sh_tlet_all_local = load_sh(shdir / "NB_Z_narrow_LET_p4.dat", scale=0.1)
    _, sh_tlet_proton_local = load_sh(shdir / "NB_Z_narrow_LET_p6.dat", scale=0.1)

    z_dose, sh_dose = load_sh(shdir / "NB_Z_narrow_dose_p2.dat", scale=1.0)

    sh_dlet_all_i = interp_to_fluka(z_sh, sh_dlet_all_local, z)
    sh_dlet_proton_i = interp_to_fluka(z_sh, sh_dlet_proton_local, z)
    sh_tlet_all_i = interp_to_fluka(z_sh, sh_tlet_all_local, z)
    sh_tlet_proton_i = interp_to_fluka(z_sh, sh_tlet_proton_local, z)

    sh_dose_i = interp_to_fluka(z_dose, sh_dose, z)
    dose_mask = sh_dose_i / np.nanmax(sh_dose_i) > 0.05
    fluka_mask = flu_total > 0.0
    mask = dose_mask & fluka_mask

    fl_t = np.where(mask, fluka_tlet, np.nan)
    fl_d = np.where(mask, fluka_dlet, np.nan)

    sh_t_all = np.where(dose_mask, sh_tlet_all_i, np.nan)
    sh_t_p = np.where(dose_mask, sh_tlet_proton_i, np.nan)
    sh_d_all = np.where(dose_mask, sh_dlet_all_i, np.nan)
    sh_d_p = np.where(dose_mask, sh_dlet_proton_i, np.nan)

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.plot(z, fl_t, linewidth=1.8, label="FLUKA scored TLET: p+d/t/He/Li")
    ax.plot(z, fl_d, linewidth=1.8, label="FLUKA scored DLET: p+d/t/He/Li")

    ax.plot(z, sh_t_all, "--", linewidth=1.5, label="SHIELD-HIT TLET local: all charged")
    ax.plot(z, sh_d_all, "--", linewidth=1.5, label="SHIELD-HIT DLET local: all charged")

    ax.plot(z, sh_t_p, ":", linewidth=1.5, label="SHIELD-HIT TLET local: protons")
    ax.plot(z, sh_d_p, ":", linewidth=1.5, label="SHIELD-HIT DLET local: protons")

    ax.set_title("{}\nCumulative LET: FLUKA scored species vs SHIELD-HIT local".format(case))
    ax.set_xlabel("z [cm]")
    ax.set_ylabel("LET [keV/um]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    ax.set_xlim(float(z[0]), float(z[-1]))
    ax.set_ylim(0.0, ylimit_from(fl_t, fl_d, sh_t_all, sh_d_all, sh_t_p, sh_d_p))

    fig.tight_layout()

    out = OUT / "{}_FLUKA_cumulative_scored_LET_vs_SHIELDHIT_local.png".format(case)
    fig.savefig(out, dpi=200)
    plt.close(fig)

    row = {
        "case": case,
        "fluence_total": float(np.sum(flu_total)),
        "LET1_total": float(np.sum(let1_total)),
        "LET2_total": float(np.sum(let2_total)),
        "global_TLET_scored_keV_um": float(np.sum(let1_total) / np.sum(flu_total)) if np.sum(flu_total) > 0.0 else "",
        "global_DLET_scored_keV_um": float(np.sum(let2_total) / np.sum(let1_total)) if np.sum(let1_total) > 0.0 else "",
        "max_depth_TLET_scored_keV_um": float(np.nanmax(fluka_tlet)),
        "max_depth_DLET_scored_keV_um": float(np.nanmax(fluka_dlet)),
    }

    for name, _, _, _ in COMPONENTS:
        row["fluence_sum_{}".format(name)] = component_flu_sums[name]
        row["LET1_sum_{}".format(name)] = component_let1_sums[name]
        row["LET2_sum_{}".format(name)] = component_let2_sums[name]

    rows.append(row)

    print(
        "{}  global_TLET={:.6g}  global_DLET={:.6g}  max_TLET={:.6g}  max_DLET={:.6g}  -> {}".format(
            case,
            row["global_TLET_scored_keV_um"],
            row["global_DLET_scored_keV_um"],
            row["max_depth_TLET_scored_keV_um"],
            row["max_depth_DLET_scored_keV_um"],
            out,
        )
    )

out_csv = Path("fluka_cumulative_scored_LET_vs_shieldhit_LOCAL_{}_summary.csv".format(TAG))

with out_csv.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print("wrote", out_csv)
print("wrote", OUT)
