# FRAGLET_COMSCW_30M milestone — 2026-06-19

## Purpose

This milestone records the current validated FLUKA vs SHIELD-HIT post-processing and plotting package for the all-8-case LET benchmark, including the new fragment scoring workflow.

The production run used the FLUKA custom fragment/LET scoring workflow with COMSCW/FLUSCW-style scoring and the corrected source/weighting setup used in the REFPLANE workflow.

## Production run

Run tag: FRAGLET_COMSCW_30M

Cases:
- plan01_field01_geoA_SOBPcent
- plan01_field01_geoB_SOBP95
- plan01_field01_geoC_SOBP74
- plan02_field01_geoD_mono
- plan03_field01_geoA_rampFull
- plan03_field02_geoA_rampFull
- plan04_field01_geoA_rampMiddle
- plan04_field02_geoA_rampMiddle

Run structure:
- 8 cases
- 10 chunks per case
- 3,000,000 primaries per chunk
- 30,000,000 primaries per case
- 240,000,000 primaries total

SLURM job:
- job id: 5422074
- job name: fluka_all8_fraglet_30M
- finished: Fri Jun 19 04:06:08 CEST 2026
- DONE lines: 80
- stderr lines: 0
- failed/error/abort scan: empty

Important run settings confirmed in the chunk inputs/logs:
- TAG=FRAGLET_COMSCW_30M
- NCHUNKS=10
- PRIMARIES=3000000
- BEAMPOS 0 0 -50 cm
- SOURCE card active
- USERWEIG active
- START 3000000

## Scored fragment species

The current cumulative FLUKA scored-fragment comparison includes:
- proton
- deuteron
- triton
- helium3
- helium4
- lithium6
- lithium7

Important caveat:

FLUKA cumulative scored LET is a selected-species cumulative LET:
p + d + t + He-3 + He-4 + Li-6 + Li-7.

SHIELD-HIT "all charged" may include additional charged fragments.
Therefore the FLUKA cumulative scored curve should not be described as a complete all-charged-particle equivalent unless all relevant charged species are explicitly scored.

## Output package

Final archive on Grendel:

/home/dewh/runs/let_benchmark/ALL_PLOTS_FOR_MAC_FRAGLET_COMSCW_30M_20260619.tar.gz

Archive size: 12M

SHA256:

11564a75939cfc25960faedf7a9988deb7fd5a571ad14dbda2c0ec34b21679cf

Verified on Mac:

~/let_benchmark_plots/ALL_PLOTS_FOR_MAC_FRAGLET_COMSCW_30M_20260619

PNG count: 88

## Plot package contents

The final package contains:
- 48 corrected normal FLUKA vs SHIELD-HIT dose/fluence/LET plots
- 24 fragment LET / DLET-composition / cumulative LET plots
- 16 fragment fluence + fragment dose plots

Subfolders:
- plots_getlet_vs_shieldhit_FRAGLET_COMSCW_30M
- plots_fluka_scored_dlet_composition_FRAGLET_COMSCW_30M
- plots_fluka_scored_dlet_composition_vs_shieldhit_LOCAL_FRAGLET_COMSCW_30M
- plots_fluka_cumulative_scored_LET_vs_shieldhit_LOCAL_FRAGLET_COMSCW_30M
- plots_fragment_fluence_dose_profiles_FRAGLET_COMSCW_30M

## Important post-processing correction

The tagged normal 48-plot builder originally summed chunked FLUKA USRBIN .lis files directly.

That was incorrect for absolute per-primary quantities, because each independent FLUKA chunk is already normalized per primary. For equal-sized chunks, chunked profiles must be averaged:

combined_per_primary = (chunk1 + chunk2 + ... + chunk10) / 10

Affected by this correction:
- absolute dose
- absolute fluence
- absolute fragment dose
- absolute fragment fluence

Not affected in practice:
- normalized dose/fluence shapes
- TLET
- DLET
- LET local/water ratios

because normalization or numerator/denominator cancellation removes the common chunk factor.

After correction, the normal 30M absolute dose metrics became consistent:

dose_abs_MAD approximately 0.0003 to 0.0006 MeV/g/primary

## Current scripts created/used

Main 30M post-processing and plotting scripts in /home/dewh/runs/let_benchmark:
- postprocess_fraglet_all8.sh
- make_getlet_doseflu_vs_shieldhit_all8_TAGGED.py
- plot_all_getlet_vs_shieldhit_TAGGED.py
- plot_fluka_scored_dlet_composition_all8.py
- plot_fluka_scored_dlet_composition_vs_shieldhit_local_all8.py
- plot_fluka_cumulative_scored_LET_vs_shieldhit_local_all8.py
- plot_fragment_fluence_dose_profiles_all8.py
- make_fraglet_all8_csvs.py
- convert_all8_fraglet_forts21_33.sh
- convert_all8_fraglet_forts41_67.sh
- slurm_all8_fraglet_10x50k.sh

## Next repository-cleanup objective

Promote the working scripts above into the maintained 2022_DCPT_LET repo copy so that the 30M workflow is reproducible from the project repository rather than only from the working run directory.

Before replacing scripts in the repo, make backups and inspect git status/git diff.
