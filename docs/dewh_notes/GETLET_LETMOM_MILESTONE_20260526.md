# GETLET LET-moment milestone — 2026-05-26

Goal:
Replace the hard-coded SHIELD-HIT/libdEdx proton stopping-power table in FLUSCW with FLUKA's internal GETLET routine.

Working test file:
- fluka_custom_sobp_letmom_titusb_getlet_test.f

Working executable:
- fluka_custom_sobp_letmom_titusb_getlet_test

Core FLUSCW logic:
```fortran
MATLET = MEDFLK(NREG,1)
EKIN   = -PLA

LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)

IF ( TITUSB(JSCRNG) .EQ. '      PHL1' ) THEN
   FLUSCW = LETW
ELSE IF ( TITUSB(JSCRNG) .EQ. '      PHL2' ) THEN
   FLUSCW = LETW * LETW
END IF
```

Important findings:
- GETLET exists in FLUKA libfluka.a as getlet_.
- The combined GETLET test file compiles with ldpmqmd.
- MEDIUM(NREG) caused a crash because IPRODC was 0 inside FLUSCW.
- MEDFLK(NREG,1) works for prompt region-to-material mapping.
- The argument order GETLET(MATLET, EKMEV, AMASS, ZERZER, IJ) compiled but produced zero PHL1/PHL2.
- The course-style argument order GETLET(IJ, EKIN, PLA, ZERZER, MATLET) produced nonzero PHL1 and PHL2.
- The hard-coded SHIELD-HIT/libdEdx EDEDX/SDEDX table has been removed from the GETLET test file.

1k mono smoke-test diagnostic:
- PHI, PHL1, and PHL2 each had 205 bins.
- Nonzero bins: PHI/PHL1/PHL2 = 184/184/184.
- Integrated TLET = PHL1/PHI  = 0.8934889014898265.
- Integrated DLET = PHL2/PHL1 = 1.9566462747254614.

Caveat:
Using MATLET = MEDFLK(NREG,1) means local-material FLUKA LET, not forced water LET.
For the mono case this means SOLWATR in SWF and PMMA in the slab.

Current status:
The GETLET-based LET-moment scorer is technically working in a 1k smoke test. It is not yet a production-statistics result.

## Forced-water LET milestone extension

Goal:
Test whether FLUKA GETLET can score proton LET using ordinary water as the stopping-power reference material, while transport still occurs in the actual geometry materials.

Input test file:
- ../fluka_plan02_field01_geoD_mono_letmom/plan02_field01_geoD_mono_letmom_titusb_names_1k_watermat_test.inp

Forced-water executable:
- fluka_custom_sobp_letmom_titusb_getlet_water_test

Water material setup:
- Added ordinary water material named WATE.
- WATE composition: 2 hydrogen + 1 oxygen.
- WATE density: 1.000 g/cm3.
- WATE was not assigned to the real phantom.
- A dummy region WATEDUM was added outside the real WORLD region.
- WATEDUM was assigned to WATE only to force FLUKA to initialize stopping-power/GETLET data for WATE.

Important finding:
- Defining WATE without assigning it to any region caused GETLET to crash with a floating-point exception.
- Assigning WATE to the dummy WATEDUM region fixed the crash.
- Therefore, for GETLET use, the water-reference material should be assigned somewhere, even if only to a dummy unreachable region.

Material numbers in the water-material mono test:
- SOLWATR = 27
- WATE    = 28
- AIR     = 29
- PMMA    = 30

Forced-water FLUSCW test logic:
- MATLET = 28
- EKIN   = -PLA
- LETW   = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)

1k mono forced-water diagnostic:
- PHI, PHL1, and PHL2 each had 205 bins.
- Nonzero bins: PHI/PHL1/PHL2 = 184/184/184.
- sum PHI  = 2.636131380639
- sum PHL1 = 2.39996703827
- sum PHL2 = 4.7896567372
- Integrated TLET_water = PHL1/PHI  = 0.9104125294727329
- Integrated DLET_water = PHL2/PHL1 = 1.995717716461886

Comparison with local-material GETLET smoke test:
- TLET_local = 0.8934889014898265
- DLET_local = 1.9566462747254614
- TLET_water = 0.9104125294727329
- DLET_water = 1.995717716461886

Interpretation:
- Local-material LET uses MATLET = MEDFLK(NREG,1), i.e. SOLWATR, PMMA, or AIR depending on the region.
- Water-reference LET uses MATLET = WATE material number, i.e. ordinary water stopping power everywhere.
- For comparison to SHIELD-HIT LET Protons in_Water, the water-reference LET is the relevant FLUKA analogue.
- For comparison to SHIELD-HIT LET Protons, the local-material LET is the relevant FLUKA analogue.

Next design goal:
Build one combined .f file and one input that score both local-material and water-reference LET moments in the same FLUKA run.

## Combined local + water LET run

Combined Fortran test file:
- fluka_custom_sobp_letmom_titusb_getlet_both_test.f

Combined executable:
- fluka_custom_sobp_letmom_titusb_getlet_both_test

Combined input:
- ../fluka_plan02_field01_geoD_mono_letmom/plan02_field01_geoD_mono_letmom_titusb_both_1k.inp

Combined USRBIN outputs:
- fort.21 = EDEP_ZN
- fort.22 = PHI_ZN
- fort.23 = PHL1_ZN = local-material proton LET moment, sum ds*LET_local
- fort.24 = PHL2_ZN = local-material proton LET squared moment, sum ds*LET_local^2
- fort.25 = PWL1_ZN = water-reference proton LET moment, sum ds*LET_water
- fort.26 = PWL2_ZN = water-reference proton LET squared moment, sum ds*LET_water^2

1k mono combined diagnostic:
- bins: 205 for all moment arrays.
- sum PHI        = 2.636131380639
- sum PHL1 local = 2.35535413147
- sum PHL2 local = 4.608594887
- sum PWL1 water = 2.39996703827
- sum PWL2 water = 4.7896567372
- TLET_local     = 0.8934889014898265
- DLET_local     = 1.9566462747254614
- TLET_water     = 0.9104125294727329
- DLET_water     = 1.995717716461886
- nonzero bins PHI/PHL1/PHL2/PWL1/PWL2 = 184/184/184/184/184

Conclusion:
The combined GETLET scorer reproduces both the local-material and water-reference LET diagnostics in one FLUKA run.

## Combined implementation audit snapshot

Fortran audit:
- File: custom_f/fluka_custom_sobp_letmom_titusb_getlet_both_test.f
- PHL1 uses MATLET = MEDFLK(NREG,1), then GETLET(IJ, EKIN, PLA, ZERZER, MATLET).
- PHL2 uses MATLET = MEDFLK(NREG,1), then GETLET(IJ, EKIN, PLA, ZERZER, MATLET), squared.
- PWL1 uses MATLET = 28, then GETLET(IJ, EKIN, PLA, ZERZER, MATLET).
- PWL2 uses MATLET = 28, then GETLET(IJ, EKIN, PLA, ZERZER, MATLET), squared.

Input audit:
- File: fluka_plan02_field01_geoD_mono_letmom/plan02_field01_geoD_mono_letmom_titusb_both_1k.inp
- WATE material is defined as ordinary water.
- WATEDUM dummy region exists outside WORLD.
- BLCKHOLE excludes WATEDUM.
- WATEDUM is assigned to WATE.
- USRBIN scorers:
  - fort.21 = EDEP_ZN
  - fort.22 = PHI_ZN
  - fort.23 = PHL1_ZN
  - fort.24 = PHL2_ZN
  - fort.25 = PWL1_ZN
  - fort.26 = PWL2_ZN
- USERWEIG is active after the USRBIN block.

Caveat:
- MATLET = 28 is valid for this mono test input because FLUKA assigned WATE material number 28.
- This should be verified or generated consistently before using the method across all benchmark cases.

## Cross-case WATE material-number check

Additional test case:
- fluka_plan01_field01_geoA_SOBPcent_letmom/plan01_field01_geoA_SOBPcent_letmom_both_test.inp

Purpose:
Check whether inserting WATE after SOLWATR gives the same FLUKA material numbering outside the mono case.

Result from FLUKA material table:
- SOLWATR = 27
- WATE    = 28
- AIR     = 29
- PMMA    = 30

Conclusion:
The WATE material-number pattern matched the mono test. This supports using MATLET = 28 for water-reference LET if the input generator consistently inserts WATE immediately after SOLWATR and assigns it to WATEDUM.

## All-case LET-moment material-block audit

Audit command:
Checked all current fluka_plan*_letmom/*_step2_prod.inp files for MATERIAL, COMPOUND, ASSIGNMA, SOLWATR, PMMA, and AIR.

Result:
All audited LET-moment production inputs have the same relevant material structure:
- SOLWATR material definition appears before AIR and PMMA assignments.
- SOLWATR compound lines immediately follow the SOLWATR MATERIAL line.
- ASSIGNMA lines assign AIR, SOLWATR, and PMMA in the same order.
- None of the current production inputs define WATE yet.

Implication:
The input generator can systematically insert WATE immediately after the SOLWATR compound block and assign it to a dummy WATEDUM region.

Expected material numbering if generated consistently:
- SOLWATR = 27
- WATE    = 28
- AIR     = 29
- PMMA    = 30

Production caveat:
The generated FLUKA .out material table should still be checked for every case before using MATLET = 28 in production LET_water comparisons.

## Generator patch validation

Patched generator:
- /home/dewh/runs/let_benchmark/make_fluka_case_from_sh12a_letmom.py

Generator now writes:
- WATEDUM dummy geometry region
- WATE ordinary-water material
- ASSIGNMA WATE WATEDUM
- PHI_ZN, PHL1_ZN, PHL2_ZN, PWL1_ZN, PWL2_ZN
- USERWEIG

Freshly generated mono test:
- fluka_plan02_field01_geoD_mono_letmom/plan02_field01_geoD_mono_letmom_generated_both_1k.inp

Result:
- The generated 1k input ran successfully with fluka_custom_sobp_letmom_titusb_getlet_both_test.
- It produced fort.21 through fort.26.
- The generated result exactly reproduced the earlier manual combined test.

Generated 1k diagnostic:
- sum PHI        = 2.636131380639
- sum PHL1 local = 2.35535413147
- sum PHL2 local = 4.608594887
- sum PWL1 water = 2.39996703827
- sum PWL2 water = 4.7896567372
- TLET_local     = 0.8934889014898265
- DLET_local     = 1.9566462747254614
- TLET_water     = 0.9104125294727329
- DLET_water     = 1.995717716461886
- nonzero bins PHI/PHL1/PHL2/PWL1/PWL2 = 184/184/184/184/184

Conclusion:
The patched LET-moment generator can now produce runnable combined local-material and water-reference proton LET inputs.

## Generator manifest cleanup

Patched generator:
- /home/dewh/runs/let_benchmark/make_fluka_case_from_sh12a_letmom.py

Cleanup performed:
- Updated the generator manifest text to match the new combined GETLET scoring scheme.
- Removed stale PHIL_ZN / PHIL2_ZN wording from the manifest text.
- Removed stale table-based Water.txt / fluscw_let wording from the manifest text.

Current manifest terminology:
- PHI_ZN  = proton track-length fluence moment, sum(ds)
- PHL1_ZN = local-material proton LET moment, sum(ds * LET_local)
- PHL2_ZN = local-material proton LET squared moment, sum(ds * LET_local^2)
- PWL1_ZN = water-reference proton LET moment, sum(ds * LET_water)
- PWL2_ZN = water-reference proton LET squared moment, sum(ds * LET_water^2)

Current reconstruction formulas:
- TLET_local = PHL1_ZN / PHI_ZN
- DLET_local = PHL2_ZN / PHL1_ZN
- TLET_water = PWL1_ZN / PHI_ZN
- DLET_water = PWL2_ZN / PWL1_ZN

Current stopping-power statement:
- LET is calculated inside the combined GETLET FLUSCW routine using FLUKA GETLET.
- Local-material LET uses MATLET = MEDFLK(NREG,1).
- Water-reference LET uses WATE, expected material number 28 when WATE is inserted immediately after SOLWATR and assigned to WATEDUM.

Status:
The LET-moment generator has been patched and validated for the mono case. Broader regeneration/testing across all benchmark cases is paused until the separate thread on additional SHIELD-HIT/FLUKA outputs is resolved.

## Merged dose/fluence + GETLET mono smoke validation

A new merged canonical Fortran source was created:

- `custom_f/fluka_custom_sobp_letmom_titusb_doseflu_getlet.f`

It was copied from:

- `fluka_custom_sobp_letmom_titusb_getlet_both_test.f`

and extended with:

- `trackr.inc`
- primary-proton filtering using `LTRACK .EQ. 1`
- primary local-material LET moments:
  - `PRI_ZN`
  - `PRL1_ZN`
  - `PRL2_ZN`
- primary water-reference LET moments:
  - `PWR1_ZN`
  - `PWR2_ZN`
- dose/fluence compatibility with:
  - `PDEP_ZN`
  - `AFLU_ZN`

The merged executable compiled successfully:

- `custom_f/fluka_custom_sobp_letmom_titusb_doseflu_getlet`

Manual mono merged smoke input:

- `fluka_plan02_field01_geoD_mono_letmom/plan02_field01_geoD_mono_letmom_doseflu_getlet_smoke.inp`

Important formatting note:

- `START          1000.0` caused FLUKA fixed-format parsing failure.
- `START       1000.0` worked.

The merged mono smoke run completed successfully and produced:

- `fort.21` = `EDEP_ZN`
- `fort.22` = `PHI_ZN`
- `fort.23` = `PHL1_ZN`
- `fort.24` = `PHL2_ZN`
- `fort.25` = `PWL1_ZN`
- `fort.26` = `PWL2_ZN`
- `fort.27` = `PDEP_ZN`
- `fort.28` = `AFLU_ZN`
- `fort.29` = `PRI_ZN`
- `fort.30` = `PRL1_ZN`
- `fort.31` = `PRL2_ZN`
- `fort.32` = `PWR1_ZN`
- `fort.33` = `PWR2_ZN`

Integrated smoke-test sums:

- `EDEP_all` = 0.0247515358814
- `PHI_protons` = 2.636131380639
- `PHL1_local` = 2.35535413147
- `PHL2_local` = 4.608594887
- `PWL1_water` = 2.39996703827
- `PWL2_water` = 4.7896567372
- `PDEP_protons` = 0.0237365357114
- `AFLU_all` = 3.0063684
- `PRI_primary` = 2.58358788
- `PRL1_primary_local` = 2.2548917
- `PRL2_primary_local` = 4.0479389
- `PWR1_primary_water` = 2.297737
- `PWR2_primary_water` = 4.2051679

Integrated ratios:

- `PDEP/EDEP` = 0.9589924368789277
- `PHI/AFLU` = 0.8768490849754143
- `PRI/PHI` = 0.9800679507004452

All-proton LET reconstruction:

- `TLET_local = PHL1_ZN / PHI_ZN = 0.8934889014898265`
- `DLET_local = PHL2_ZN / PHL1_ZN = 1.9566462747254614`
- `TLET_water = PWL1_ZN / PHI_ZN = 0.9104125294727329`
- `DLET_water = PWL2_ZN / PWL1_ZN = 1.995717716461886`

Primary-proton LET reconstruction:

- `TLET_primary_local = PRL1_ZN / PRI_ZN = 0.8727753050149778`
- `DLET_primary_local = PRL2_ZN / PRL1_ZN = 1.7951810723326536`
- `TLET_primary_water = PWR1_ZN / PRI_ZN = 0.8893589483784078`
- `DLET_primary_water = PWR2_ZN / PWR1_ZN = 1.830134562832909`

Conclusion:

The merged `doseflu_getlet` source is now validated for mono 1k smoke statistics. It combines GETLET local/water LET, primary-proton filtering, proton dose filtering, all-particle fluence, and all relevant LET moment families in one executable. The next step is to encode this full scorer layout into the input generator.

## Generator validation for merged doseflu_getlet layout

The patched generator:

- `/home/dewh/runs/let_benchmark/make_fluka_case_from_sh12a_letmom.py`

was updated to produce the full merged scorer layout:

- `fort.21` = `EDEP_ZN`
- `fort.22` = `PHI_ZN`
- `fort.23` = `PHL1_ZN`
- `fort.24` = `PHL2_ZN`
- `fort.25` = `PWL1_ZN`
- `fort.26` = `PWL2_ZN`
- `fort.27` = `PDEP_ZN`
- `fort.28` = `AFLU_ZN`
- `fort.29` = `PRI_ZN`
- `fort.30` = `PRL1_ZN`
- `fort.31` = `PRL2_ZN`
- `fort.32` = `PWR1_ZN`
- `fort.33` = `PWR2_ZN`

It also preserves:

- `WATEDUM`
- `WATE`
- `ASSIGNMA WATE WATEDUM`
- `AUXSCORE USRBIN PROTON PDEP`
- `USERWEIG`

The smoke `START` formatting was corrected from:

- `START          10.0`

to:

- `START       10.0`

because the former caused FLUKA fixed-format parsing problems when manually testing `START          1000.0`.

The mono case was regenerated with:

```bash
python3 make_fluka_case_from_sh12a_letmom.py plan02_field01_geoD_mono