# FLUKA CERN data

This directory contains FLUKA input and/or result material for the DCPT LET benchmark.

The spot lists can be read by compiling the FLUKA SOBP source sampler from:

https://github.com/DataMedSci/au-fluka-tools/tree/master/fluka_sobp_source

and linking it to FLUKA, then specifying the corresponding SOURCE card in the FLUKA input file.

## Important dose-scorer interpretation note

The difference between all-particle dose and proton-filtered dose must not be interpreted as charged-fragment dose.

In the FLUKA scoring used here, the relevant quantities are:

- DOSE_ZN: all-particle dose
- PDOSE_ZN: proton-filtered dose
- P1DO_ZN or PRDO_ZN: primary-proton-filtered dose, where available
- explicit fragment-dose scorers: deuteron, triton, helium-3, helium-4, lithium-6, lithium-7, and generic HEAVYION, where available

A diagnostic decomposition of the FLUKA all-minus-proton dose residual, DOSE_ZN - PDOSE_ZN, gave the following category attribution:

| Component | Fraction of DOSE_ZN - PDOSE_ZN |
|---|---:|
| D/T/He-3/He-4 dose | 20.207% |
| Li-6 + Li-7 dose | 0.000% |
| generic transported HEAVYION dose | 1.492% |
| neutron dose | 0.448% |
| photon dose | 0.002% |
| electron dose | 77.849% |
| positron dose | 0.000% |
| remaining residual | 0.001% |

Therefore, DOSE_ZN - PDOSE_ZN is a broad non-proton-attributed dose residual, dominated in this diagnostic by electron-attributed dose. It is not a charged-fragment dose estimator.

Fragment-dose analysis should use the explicit fragment-dose scorers, not the difference between all-particle dose and proton-filtered dose.
