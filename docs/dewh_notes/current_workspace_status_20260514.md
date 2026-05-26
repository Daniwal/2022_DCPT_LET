# DEWH FLUKA benchmark workspace status - 2026-05-14

Repository branch:
dewh/fluka-let-benchmark-workspace

Remotes:
origin   = https://github.com/Daniwal/2022_DCPT_LET.git
upstream = https://github.com/APTG/2022_DCPT_LET.git

Upstream push has been disabled intentionally.

Current purpose:
This workspace organizes MSc thesis FLUKA vs SHIELD-HIT benchmark work, including FLUKA input cards, SHIELD-HIT-derived FLUKA cases, source sampler files, LET-moment input variants, fragmentation test inputs, workflow notes, and helper scripts.

Production outputs and binaries are intentionally excluded from Git.

FLAIR status:
FLAIR is installed on the Mac and can open copied FLUKA input files for card inspection.
The macOS geoviewer plugin was compiled and imported in a small test, but crashes FLAIR during normal startup, so it is disabled for now.
FLUKA execution remains on Grendel.

LET-moment status:
PHI_ZN   = sum(ds)
PHIL_ZN  = sum(ds * LET)
PHIL2_ZN = sum(ds * LET^2)

LET_T = PHIL_ZN / PHI_ZN
LET_D = PHIL2_ZN / PHIL_ZN

A combined custom Fortran file exists outside the repository working tree:
/home/dewh/runs/let_benchmark/custom_f/fluka_custom_sobp_letmom.f

It combines SOBP source sampling and FLUSCW LET-moment scoring.

Git hygiene check:
A check was run for suspicious Git-visible outputs/binaries.
Result: No suspicious files found.

Not ready to commit yet:
1. Review final file set.
2. Separate fragmentation files from baseline LET/dose files.
3. Decide combined custom .f strategy.
4. Check generated input folders for completeness.
5. Review .gitignore one final time.

Additional FLAIR limitation:
FLAIR on Mac can open and inspect input cards, but Save/Save As froze with a _tkinter.TclError and did not create the test file. Therefore FLAIR should currently be treated as inspection-only. All actual edits should be done in VS Code / Grendel, not saved from FLAIR.
