#!/usr/bin/env bash
#SBATCH --job-name=fluka_all8_fraglet_10x50k
#SBATCH --partition=q28
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --output=slurm_all8_fraglet_10x50k_%j.out
#SBATCH --error=slurm_all8_fraglet_10x50k_%j.err

set -euo pipefail

module purge
module load gcc/9.3.0

export FLUPRO=/home/dewh/software/fluka_deb/usr/local/fluka
export FLUKADATA="${FLUPRO}/data"

ROOT=/home/dewh/runs/let_benchmark

CASES=(
plan01_field01_geoA_SOBPcent
plan01_field01_geoB_SOBP95
plan01_field01_geoC_SOBP74
plan02_field01_geoD_mono
plan03_field01_geoA_rampFull
plan03_field02_geoA_rampFull
plan04_field01_geoA_rampMiddle
plan04_field02_geoA_rampMiddle
)

NCHUNKS="${NCHUNKS:-10}"
PRIMARIES="${PRIMARIES:-50000}"
MAX_PARALLEL="${MAX_PARALLEL:-20}"
TAG="${TAG:-FRAGLET_COMSCW_10x50k}"

echo "==== job started at $(date) on $(hostname) ===="
echo "ROOT=${ROOT}"
echo "NCHUNKS=${NCHUNKS}"
echo "PRIMARIES=${PRIMARIES}"
echo "MAX_PARALLEL=${MAX_PARALLEL}"
echo "TAG=${TAG}"

running_jobs=0

run_chunk () {
local case_name="$1"
local task="$2"

local cdir="${ROOT}/fluka_${case_name}_letmom"
local base="${case_name}_letmom_step2_prod.inp"
local inp="${case_name}_letmom_step2_prod_${TAG}_chunk${task}_run.inp"
local exe="./fluka_custom_sobp_letmom_titusb_doseflu_fraglet"

local seed
seed=$((3000000 + 1000 * task + $(echo "${case_name}" | cksum | awk '{print $1 % 900}')))

cd "${cdir}"

cp "${base}" "${inp}"

python3 - <<PY
from pathlib import Path
import re

p = Path("${inp}")
s = p.read_text()

seed = float(${seed})
primaries = float(${PRIMARIES})

s = re.sub(
r"^RANDOMIZ\s+1.0\s+[-+0-9.Ee]+.*$",
"RANDOMIZ         1.0  {seed:10.1f}".format(seed=seed),
s,
flags=re.MULTILINE,
)

s = re.sub(
r"^START\s+[-+0-9.Ee]+.*$",
"START     {primaries:10.1f}".format(primaries=primaries),
s,
flags=re.MULTILINE,
)

p.write_text(s)
PY

echo "START case=${case_name} chunk=${task} seed=${seed} primaries=${PRIMARIES}"
grep -nE "USERWEIG|RANDOMIZ|START" "${inp}"

"${FLUPRO}/bin/rfluka" -N 0 -M 1 -e "${exe}" "${inp%.inp}"

out="${case_name}_letmom_step2_prod_${TAG}_chunk${task}_run001.out"

if grep -q "End of FLUKA v4 run" "${out}"; then
echo "DONE case=${case_name} chunk=${task}"
else
echo "FAILED_OR_INCOMPLETE case=${case_name} chunk=${task}"
return 1
fi
}

export -f run_chunk
export ROOT FLUPRO FLUKADATA TAG PRIMARIES

for case_name in "${CASES[@]}"; do
for task in $(seq 1 "${NCHUNKS}"); do
run_chunk "${case_name}" "${task}" &

running_jobs=$((running_jobs + 1))

if [[ "${running_jobs}" -ge "${MAX_PARALLEL}" ]]; then
  wait -n
  running_jobs=$((running_jobs - 1))
fi

done
done

wait

echo "==== all chunks finished at $(date) ===="
