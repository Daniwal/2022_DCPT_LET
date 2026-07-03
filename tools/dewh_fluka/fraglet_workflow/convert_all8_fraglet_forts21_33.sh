#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-/home/dewh/runs/let_benchmark}"
TAG="${TAG:-FRAGLET_COMSCW_10x50k}"
NCHUNKS="${NCHUNKS:-1}"

cases=(
plan01_field01_geoA_SOBPcent
plan01_field01_geoB_SOBP95
plan01_field01_geoC_SOBP74
plan02_field01_geoD_mono
plan03_field01_geoA_rampFull
plan03_field02_geoA_rampFull
plan04_field01_geoA_rampMiddle
plan04_field02_geoA_rampMiddle
)

forts=(21 22 23 24 25 26 27 28 29 30 31 32 33 34)

for case in "${cases[@]}"; do
d="${ROOT}/fluka_${case}_letmom"
cd "${d}"

echo
echo "===== CONVERT ${case} ====="

for chunk in $(seq 1 "${NCHUNKS}"); do
base="${case}_letmom_step2_prod_${TAG}_chunk${chunk}_run001"
echo "---- chunk ${chunk}: ${base} ----"

for fort in "${forts[@]}"; do
  in_file="${base}_fort.${fort}"
  out_file="${case}_${TAG}_chunk${chunk}_fort${fort}.lis"

  if [[ ! -f "${in_file}" ]]; then
    echo "MISSING ${in_file}"
    continue
  fi

  printf "%s\n%s\n" "${in_file}" "${out_file}" | "${FLUPRO}/bin/usbrea" >/dev/null
  echo "converted fort.${fort} -> ${out_file}"
done
done
done
