# Custom FLUKA Fortran consolidation milestone — 2026-05-14

Goal:
Collect the custom FLUKA work into one master Fortran source file instead of relying on scattered separate files.

Canonical source files:
- custom_f/source_sampler.f
- custom_f/fluscw_let.f

Canonical source hashes:
- source_sampler.f: e00e1f41414ac021a0898d1018832ba995d65a5d0e0b4273a8046799bf0039f1
- fluscw_let.f: 633f2e727cddafbeb118e26e2a76c94c3735b9fecb9c15df30c9ff7cad613a49

Combined source:
- custom_f/fluka_custom_sobp_letmom.f

Combined source hash:
- c4adc3ea524f93f414523f667c958a682c93913a452368f4fb56d41659fe1202

Compiled executable:
- custom_f/fluka_custom_sobp_letmom

Build command:
$FLUPRO/bin/ldpmqmd -o fluka_custom_sobp_letmom fluka_custom_sobp_letmom.f

Validation smoke tests:
1. plan02_field01_geoD_mono_letmom
   - Input: plan02_field01_geoD_mono_letmom_combined_smoke.inp
   - Result: End of FLUKA run
   - SOBP SOURCE READING ../sobpcln
   - SOBP NUMBER OF COLUMNS 7
   - SOBP SOURCE beamlets found: 323
   - SOBP POINT-LIKE VIRTUAL SOURCE
   - fort.21–24 produced

2. plan03_field01_geoA_rampFull_letmom
   - Input: plan03_field01_geoA_rampFull_letmom_combined_smoke.inp
   - Result: End of FLUKA run
   - SOBP SOURCE READING ../sobpcln
   - SOBP NUMBER OF COLUMNS 7
   - SOBP SOURCE beamlets found: 9218
   - SOBP POINT-LIKE VIRTUAL SOURCE
   - fort.21–24 produced

Conclusion:
The SOBP source sampler and LET-moment FLUSCW routine can be combined into one master Fortran file and compiled into one executable without linker/runtime problems.

## All-8 smoke validation with central executable

After copying the same executable into all 8 `_letmom` case folders, the script

- run_all_custom_f_smokes.sh

was used to run all combined smoke inputs with:

$FLUPRO/bin/rfluka -N 0 -M 1 -e ./fluka_custom_sobp_letmom <case>_letmom_combined_smoke.inp

All 8 cases completed with:
- End of FLUKA run
- fort.21–24 produced
- SOBP SOURCE READING ../sobpcln
- SOBP NUMBER OF COLUMNS 7
- SOBP POINT-LIKE VIRTUAL SOURCE

Validated beamlet counts:
- plan01_field01_geoA_SOBPcent: 6069
- plan01_field01_geoB_SOBP95: 6069
- plan01_field01_geoC_SOBP74: 6069
- plan02_field01_geoD_mono: 323
- plan03_field01_geoA_rampFull: 9218
- plan03_field02_geoA_rampFull: 9218
- plan04_field01_geoA_rampMiddle: 7708
- plan04_field02_geoA_rampMiddle: 9218

Log file:
- custom_f_smoke_all8_20260514.log

Conclusion:
The one-file custom FLUKA source/executable is validated for all 8 benchmark cases at smoke-test level.

## TITUSB / scorer-title investigation

Question:
Can the LET-moment FLUSCW logic identify PHI/PHIL/PHIL2 by scorer title instead of by fragile JSCRNG order?

Findings:
- `scohlp.inc` provides ISCRNG and JSCRNG.
- `usrbin.inc` provides USRBIN title array:
  - TITUSB(MXUSBN)
- Inside FLUSCW, TITUSB(JSCRNG) is accessible at runtime for USRBIN track-length scoring.

Important limitation:
The original scorer names were:
- PHI_ZN
- PHIL_ZN
- PHIL2_ZN

Runtime TITUSB values showed only the first 4 effective title characters:
- PHI_ZN   -> "PHI_"
- PHIL_ZN  -> "PHIL"
- PHIL2_ZN -> "PHIL"

Therefore, the original PHIL_ZN and PHIL2_ZN names are not distinguishable by TITUSB.

Robust naming solution:
Use unique first-four-character scorer names:
- PHI_ZN
- PHL1_ZN
- PHL2_ZN

Runtime TITUSB values then become:
- PHI_ZN   -> "PHI_"
- PHL1_ZN  -> "PHL1"
- PHL2_ZN  -> "PHL2"

Clean title-based source:
- custom_f/fluka_custom_sobp_letmom_titusb.f

Hashes:
- fluka_custom_sobp_letmom_titusb.f: f1500a5ccb8ee58152c3959142285529e3d91b15731107cb5355621817c5cfda
- fluka_custom_sobp_letmom_titusb:   29925bb741ac2a69d61977f8e3833d09bfd7b73dba2ec28580387707f308f105

Title-based FLUSCW logic:
- If TITUSB(JSCRNG) == "      PHL1", apply LET weight.
- If TITUSB(JSCRNG) == "      PHL2", apply LET^2 weight.
- Otherwise leave FLUSCW = 1.

Validation:
A 1000-primary mono test using PHL1_ZN/PHL2_ZN completed successfully and produced fort.21-24.
Converted moment sums:
- PHI  sum  = 2.6361313806389988
- PHL1 sum  = 2.4088830322999986
- PHL2 sum  = 4.64034375
- PHL1/PHI  = 0.9137947562067581
- PHL2/PHL1 = 1.9263466460508902

Interpretation:
These low-statistics integrated ratios are not physics results, but they confirm that PHL1 and PHL2 are being weighted differently by the title-based FLUSCW routine.

Fragmentation / other scorer relevance:
TITUSB is specific to USRBIN scorers. Other scoring families have their own title arrays:
- USRBIN    -> TITUSB  from usrbin.inc
- USRTRACK  -> TITUTC  from usrtrc.inc
- USRBDX    -> TITUSX  from usrbdx.inc
- USRYIELD  -> TITUYL  from usryld.inc
- RESNUCLEi -> TIURSN  from usrsnc.inc

Conclusion:
Title-array based scorer identification is useful, but must be matched to the scoring family. For LET-moment USRBIN scoring, TITUSB with unique first-four-character names is more robust than hard-coding JSCRNG=3/4.

## All-8 TITUSB smoke validation

The robust title-based LET-moment workflow was propagated to all 8 `_letmom` benchmark cases.

New input naming:
- <case>_letmom_titusb_smoke.inp
- <case>_letmom_titusb_prod.inp

Scorer titles now use unique first-four-character names:
- PHI_ZN
- PHL1_ZN
- PHL2_ZN

Executable:
- fluka_custom_sobp_letmom_titusb

Executable SHA256:
- 29925bb741ac2a69d61977f8e3833d09bfd7b73dba2ec28580387707f308f105

All-8 smoke command:
./run_all_titusb_smokes.sh | tee titusb_smoke_all8_20260514.log

All 8 cases completed with:
- End of FLUKA run
- SOBP SOURCE READING ../sobpcln
- SOBP NUMBER OF COLUMNS 7
- SOBP POINT-LIKE VIRTUAL SOURCE
- fort.21–24 produced

Validated beamlet counts:
- plan01_field01_geoA_SOBPcent: 6069
- plan01_field01_geoB_SOBP95: 6069
- plan01_field01_geoC_SOBP74: 6069
- plan02_field01_geoD_mono: 323
- plan03_field01_geoA_rampFull: 9218
- plan03_field02_geoA_rampFull: 9218
- plan04_field01_geoA_rampMiddle: 7708
- plan04_field02_geoA_rampMiddle: 9218

Conclusion:
The title-based TITUSB workflow is smoke-validated across all 8 benchmark cases and is more robust than hard-coding JSCRNG=3/4.


## TITUSB production validation against previous JSCRNG workflow

All 8 benchmark cases were rerun with the title-based TITUSB workflow using:

- `<case>_letmom_titusb_prod.inp`
- `fluka_custom_sobp_letmom_titusb`
- `START 100000.0`

Each case produced:

- `titusb_prod_fort21.lis` → EDEP_ZN
- `titusb_prod_fort22.lis` → PHI_ZN
- `titusb_prod_fort23.lis` → PHL1_ZN
- `titusb_prod_fort24.lis` → PHL2_ZN

The TITUSB `.lis` files were converted to:

- `<case>_letmom_titusb_prod_masked.csv`

A direct numerical comparison against the previous validated JSCRNG/order-based CSVs gave:

- max_abs_csv_difference = 0 for all 8 cases

Conclusion:
The title-based TITUSB workflow reproduces the previous validated LET-moment production results exactly, while avoiding fragile dependence on JSCRNG scorer order. This should be treated as the preferred thesis workflow going forward.