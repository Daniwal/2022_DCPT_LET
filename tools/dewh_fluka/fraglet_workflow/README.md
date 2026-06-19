# FLUKA fragment LET workflow

This folder contains the maintained/main FLUKA post-processing and plotting workflow for the all-8-case LET benchmark with fragment scoring.

The dated folder

    tools/dewh_fluka/fraglet_30M_workflow_20260619

is the frozen milestone snapshot that produced the validated 30M plot package.

This folder

    tools/dewh_fluka/fraglet_workflow

is the stable workflow folder to keep developing and reusing.

## Validated milestone

Validated run tag:

    FRAGLET_COMSCW_30M

Run structure:

    8 benchmark cases
    10 chunks per case
    3,000,000 primaries per chunk
    30,000,000 primaries per case
    240,000,000 total primaries

Final validated archive:

    /home/dewh/runs/let_benchmark/ALL_PLOTS_FOR_MAC_FRAGLET_COMSCW_30M_20260619.tar.gz

SHA256:

    11564a75939cfc25960faedf7a9988deb7fd5a571ad14dbda2c0ec34b21679cf

Final Mac folder:

    ~/let_benchmark_plots/ALL_PLOTS_FOR_MAC_FRAGLET_COMSCW_30M_20260619

Final plot count:

    88 PNGs

## Important scientific caveat

The cumulative FLUKA scored-fragment LET currently includes selected charged species:

    proton
    deuteron
    triton
    helium3
    helium4
    lithium6
    lithium7

Therefore the FLUKA cumulative scored-fragment curve is a selected-species cumulative LET curve. It should not automatically be described as a full all-charged-particle equivalent to SHIELD-HIT all charged unless all relevant charged species are explicitly scored.

## Configurable paths

The scripts keep Grendel defaults but allow path overrides.

Environment variables:

    ROOT    FLUKA working directory containing fluka_<case>_letmom folders
    SHROOT  SHIELD-HIT result directory containing one folder per case
    FLUPRO  FLUKA installation path
    TAG     run tag, e.g. FRAGLET_COMSCW_30M
    NCHUNKS number of chunks per case

Typical Grendel defaults:

    ROOT=/home/dewh/runs/let_benchmark
    SHROOT=/home/dewh/2022_DCPT_LET/data/sh12a/results
    FLUPRO=/home/dewh/software/fluka_deb/usr/local/fluka

## Typical post-processing command

    cd /home/dewh/runs/let_benchmark
    TAG=FRAGLET_COMSCW_30M NCHUNKS=10 bash /home/dewh/runs/let_benchmark/github_work/2022_DCPT_LET/tools/dewh_fluka/fraglet_workflow/postprocess_fraglet_all8.sh

## Chunk averaging rule

FLUKA USRBIN output from each independent chunk is already normalized per primary. Equal-primary chunks must therefore be averaged, not summed, when constructing absolute per-primary dose or fluence quantities.

Correct rule:

    combined_per_primary = (chunk1 + chunk2 + ... + chunkN) / N

This matters for:

    absolute dose
    absolute fluence
    fragment dose
    fragment fluence

LET ratios such as TLET and DLET are mostly unaffected by the common chunk factor because the factor cancels between numerator and denominator.

## Main scripts

    postprocess_fraglet_all8.sh
    convert_all8_fraglet_forts21_33.sh
    convert_all8_fraglet_forts41_67.sh
    make_getlet_doseflu_vs_shieldhit_all8_TAGGED.py
    plot_all_getlet_vs_shieldhit_TAGGED.py
    make_fraglet_all8_csvs.py
    plot_fluka_scored_dlet_composition_all8.py
    plot_fluka_scored_dlet_composition_vs_shieldhit_local_all8.py
    plot_fluka_cumulative_scored_LET_vs_shieldhit_local_all8.py
    plot_fragment_fluence_dose_profiles_all8.py
    slurm_all8_fraglet_10x50k.sh
