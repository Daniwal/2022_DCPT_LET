#!/usr/bin/env bash
set -euo pipefail

module purge
module load gcc/9.3.0
module load anaconda3/2021.05

export FLUPRO=/home/dewh/software/fluka_deb/usr/local/fluka
export FLUKADATA="${FLUPRO}/data"

export TAG="${TAG:-FRAGLET_COMSCW_1M}"
export NCHUNKS="${NCHUNKS:-10}"

echo "TAG=${TAG}"
echo "NCHUNKS=${NCHUNKS}"
echo

echo "=== Convert proton / primary LET forts 21-33 ==="
bash convert_all8_fraglet_forts21_33.sh

echo
echo "=== Convert fragment forts 41-67 ==="
bash convert_all8_fraglet_forts41_67.sh

echo
echo "=== Build fragment CSVs ==="
python3 make_fraglet_all8_csvs.py

echo
echo "=== Plot FLUKA scored DLET composition ==="
python3 plot_fluka_scored_dlet_composition_all8.py

echo
echo "=== Plot FLUKA scored DLET composition vs SHIELD-HIT local ==="
python3 plot_fluka_scored_dlet_composition_vs_shieldhit_local_all8.py

echo
echo "=== Plot cumulative FLUKA scored LET vs SHIELD-HIT local ==="
python3 plot_fluka_cumulative_scored_LET_vs_shieldhit_local_all8.py

echo
echo "=== Done ==="
