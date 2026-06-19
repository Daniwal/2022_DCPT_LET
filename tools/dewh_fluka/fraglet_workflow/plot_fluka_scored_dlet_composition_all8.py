#!/usr/bin/env python3

import csv
import os
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


TAG = os.environ.get("TAG", "FRAGLET_COMSCW_10x50k")
NCHUNKS = int(os.environ.get("NCHUNKS", "1"))
NZ = 205

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

COMPONENTS = [
    ("proton",   23, 24),
    ("deuteron", 52, 56),
    ("triton",   53, 57),
    ("helium3",  54, 58),
    ("helium4",  55, 59),
    ("lithium6", 50, 60),
    ("lithium7", 51, 61),
]


def read_usrbin_lis_values(path, nz=NZ):
    text = Path(path).read_text(errors="replace")
    lines = text.splitlines()

    start = None
    for i, line in enumerate(lines):
        if "Data follow" in line:
            start = i + 1
            break

    if start is None:
        raise RuntimeError("Could not find 'Data follow' in {}".format(path))

    values = []
    for line in lines[start:]:
        if "Percentage errors follow" in line:
            break
        for token in line.split():
            try:
                values.append(float(token.replace("D", "E")))
            except ValueError:
                pass

    if len(values) < nz:
        raise RuntimeError("{} has only {} values, expected {}".format(path, len(values), nz))

    return np.array(values[:nz], dtype=float)


def z_axis_for_case(case):
    if case == "plan02_field01_geoD_mono":
        zmin = -2.25
        zmax = 18.25
    else:
        zmin = -10.25
        zmax = 10.25

    dz = (zmax - zmin) / NZ
    return zmin + (np.arange(NZ) + 0.5) * dz


def safe_divide(num, den):
    out = np.full_like(num, np.nan, dtype=float)
    mask = den > 0.0
    out[mask] = num[mask] / den[mask]
    return out


def main():
    outdir = Path("plots_fluka_scored_dlet_composition_{}".format(TAG))
    outdir.mkdir(exist_ok=True)

    csv_path = Path("fluka_scored_dlet_composition_{}_all8.csv".format(TAG))
    rows = []

    for case in CASES:
        cdir = Path("fluka_{}_letmom".format(case))
        z = z_axis_for_case(case)

        let1_by_component = {}
        let2_by_component = {}

        for name, fort_let1, fort_let2 in COMPONENTS:
            let1_sum = np.zeros(NZ)
            let2_sum = np.zeros(NZ)

            for chunk in range(1, NCHUNKS + 1):
                f1 = cdir / "{}_{}_chunk{}_fort{}.lis".format(case, TAG, chunk, fort_let1)
                f2 = cdir / "{}_{}_chunk{}_fort{}.lis".format(case, TAG, chunk, fort_let2)

                let1_sum += read_usrbin_lis_values(f1)
                let2_sum += read_usrbin_lis_values(f2)

            let1_by_component[name] = let1_sum
            let2_by_component[name] = let2_sum

        denom = np.zeros(NZ)
        total_let2 = np.zeros(NZ)

        for name, _, _ in COMPONENTS:
            denom += let1_by_component[name]
            total_let2 += let2_by_component[name]

        total_dlet = safe_divide(total_let2, denom)

        component_curves = []
        component_labels = []

        for name, _, _ in COMPONENTS:
            comp = safe_divide(let2_by_component[name], denom)
            component_curves.append(np.nan_to_num(comp, nan=0.0))
            component_labels.append(name)

        for i in range(NZ):
            row = {
                "case": case,
                "z_bin": i,
                "z_cm": z[i],
                "denominator_LET1_sum": denom[i],
                "scored_DLET_keV_per_um": "" if np.isnan(total_dlet[i]) else total_dlet[i],
            }

            for j, name in enumerate(component_labels):
                val = component_curves[j][i]
                row["component_{}_keV_per_um".format(name)] = val

            rows.append(row)

        fig, ax = plt.subplots(figsize=(9, 5))

        ax.stackplot(z, component_curves, labels=component_labels, alpha=0.85)
        ax.plot(z, np.nan_to_num(total_dlet, nan=0.0), linewidth=1.5, label="sum")

        ax.set_title("{}\nFLUKA scored DLET composition: p + d/t/He/Li".format(case))
        ax.set_xlabel("z [cm]")
        ax.set_ylabel("DLET contribution [keV/um]")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=8, ncol=2)
        ax.set_xlim(float(z[0]), float(z[-1]))

        fig.tight_layout()

        out_png = outdir / "{}_FLUKA_scored_DLET_composition.png".format(case)
        fig.savefig(out_png, dpi=200)
        plt.close(fig)

        finite = total_dlet[np.isfinite(total_dlet)]
        max_dlet = float(np.max(finite)) if finite.size else 0.0
        print("{}  max_scored_DLET={:.6g}  -> {}".format(case, max_dlet, out_png))

    if rows:
        fieldnames = list(rows[0].keys())
        with csv_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    print("wrote", csv_path)


if __name__ == "__main__":
    main()
