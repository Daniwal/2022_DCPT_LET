#!/usr/bin/env python3

from pathlib import Path
import os
import csv
import math

import numpy as np


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

SPECIES = {
    "deuteron": {"flu": 41, "let1": 52, "let2": 56, "dose": 62},
    "triton": {"flu": 42, "let1": 53, "let2": 57, "dose": 63},
    "helium3": {"flu": 43, "let1": 54, "let2": 58, "dose": 64},
    "helium4": {"flu": 46, "let1": 55, "let2": 59, "dose": 65},
    "lithium6": {"flu": 48, "let1": 50, "let2": 60, "dose": 66},
    "lithium7": {"flu": 49, "let1": 51, "let2": 61, "dose": 67},
}


def read_usrbin_lis_values(path: Path, nz: int = NZ) -> np.ndarray:
    text = path.read_text(errors="ignore")
    lines = text.splitlines()

    start = None
    for i, line in enumerate(lines):
        if "Data follow" in line:
            start = i + 1
            break

    if start is None:
        raise RuntimeError(f"No 'Data follow' marker found in {path}")

    vals = []
    for line in lines[start:]:
        if "Percentage errors follow" in line:
            break

        parts = line.split()
        if not parts:
            continue

        try:
            vals.extend(float(x.replace("D", "E")) for x in parts)
        except ValueError:
            continue

    arr = np.array(vals, dtype=float)

    if arr.size < nz:
        raise RuntimeError(
            f"Too few data bins in {path}: found {arr.size}, expected {nz}"
        )

    return arr[:nz]


def safe_divide(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    out = np.full_like(numerator, np.nan, dtype=float)
    mask = denominator > 0.0
    out[mask] = numerator[mask] / denominator[mask]
    return out


def infer_z_axis(case: str) -> np.ndarray:
    zmin = -2.25
    zmax = 18.25 if case == "plan02_field01_geoD_mono" else 10.25
    dz = (zmax - zmin) / NZ
    return zmin + (np.arange(NZ) + 0.5) * dz


def finite_or_blank(x: float):
    if isinstance(x, float) and not math.isfinite(x):
        return ""
    return x
def build_rows():
    summary_rows = []
    profile_rows = []

    print("Building fragment LET/dose CSV rows")
    print("=" * 100)
    print(f"TAG = {TAG}")
    print(f"NCHUNKS = {NCHUNKS}")

    for case in CASES:
        print()
        print(case)

        case_dir = Path(f"fluka_{case}_letmom")
        z_cm = infer_z_axis(case)

        for species_name, mapping in SPECIES.items():
            flu = np.zeros(NZ, dtype=float)
            let1 = np.zeros(NZ, dtype=float)
            let2 = np.zeros(NZ, dtype=float)
            dose = np.zeros(NZ, dtype=float)

            for chunk in range(1, NCHUNKS + 1):
                flu_path = case_dir / f"{case}_{TAG}_chunk{chunk}_fort{mapping['flu']}.lis"
                let1_path = case_dir / f"{case}_{TAG}_chunk{chunk}_fort{mapping['let1']}.lis"
                let2_path = case_dir / f"{case}_{TAG}_chunk{chunk}_fort{mapping['let2']}.lis"
                dose_path = case_dir / f"{case}_{TAG}_chunk{chunk}_fort{mapping['dose']}.lis"

                flu += read_usrbin_lis_values(flu_path)
                let1 += read_usrbin_lis_values(let1_path)
                let2 += read_usrbin_lis_values(let2_path)
                dose += read_usrbin_lis_values(dose_path)

            tlet = safe_divide(let1, flu)
            dlet = safe_divide(let2, let1)

            flu_sum = float(np.sum(flu))
            let1_sum = float(np.sum(let1))
            let2_sum = float(np.sum(let2))
            dose_sum = float(np.sum(dose))

            tlet_global = let1_sum / flu_sum if flu_sum > 0.0 else float("nan")
            dlet_global = let2_sum / let1_sum if let1_sum > 0.0 else float("nan")

            summary_rows.append(
                {
                    "case": case,
                    "species": species_name,
                    "fluence_sum_cm2_per_primary": flu_sum,
                    "dose_sum_GeV_per_g_per_primary": dose_sum,
                    "LET1_sum": let1_sum,
                    "LET2_sum": let2_sum,
                    "TLET_global_keV_per_um": finite_or_blank(tlet_global),
                    "DLET_global_keV_per_um": finite_or_blank(dlet_global),
                    "nonzero_fluence_bins": int(np.count_nonzero(flu)),
                    "nonzero_dose_bins": int(np.count_nonzero(dose)),
                }
            )

            print(
                f"  {species_name:8s} "
                f"flu={flu_sum:12.5e} "
                f"dose={dose_sum:12.5e} "
                f"TLET={tlet_global:9.3g} "
                f"DLET={dlet_global:9.3g}"
            )

            for i in range(NZ):
                profile_rows.append(
                    {
                        "case": case,
                        "species": species_name,
                        "z_bin": i + 1,
                        "z_cm": float(z_cm[i]),
                        "fluence_cm2_per_primary": float(flu[i]),
                        "dose_GeV_per_g_per_primary": float(dose[i]),
                        "LET1": float(let1[i]),
                        "LET2": float(let2[i]),
                        "TLET_keV_per_um": finite_or_blank(float(tlet[i])),
                        "DLET_keV_per_um": finite_or_blank(float(dlet[i])),
                    }
                )

    return summary_rows, profile_rows
def write_csv(path: Path, rows: list) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write for {path}")

    fieldnames = list(rows[0].keys())

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    summary_rows, profile_rows = build_rows()

    summary_path = Path(f"fraglet_{TAG}_all8_summary.csv")
    profile_path = Path(f"fraglet_{TAG}_all8_depth_profiles.csv")

    write_csv(summary_path, summary_rows)
    write_csv(profile_path, profile_rows)

    print()
    print("Wrote:")
    print(f"  {summary_path}")
    print(f"  {profile_path}")


if __name__ == "__main__":
    main()