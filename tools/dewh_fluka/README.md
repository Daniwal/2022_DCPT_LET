# dewh_fluka tools

This folder contains local FLUKA/SHIELD-HIT workflow tools developed for the MSc thesis benchmark.

## Current status

`make_fluka_case_from_sh12a.py` is currently a Grendel-local prototype.

It is useful for documenting and reproducing the current workflow, but it is not yet fully portable because it contains hard-coded paths such as:

- `/home/dewh/2022_DCPT_LET`
- `/home/dewh/runs/let_benchmark`
- `/home/dewh/au-fluka-tools`

## Before treating as a general repository tool

The script should be refactored to accept command-line arguments for:

- SHIELD-HIT input root
- FLUKA output root
- FLUKA template input
- source sampler executable
- source sampler Fortran source

Until then, use it as a documented local thesis workflow script, not as a polished public tool.
