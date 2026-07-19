# VODP M0/M1 PySCF Workflow

This contains the minimal files needed to run the PySCF part of the VODP excited state project. The M0 and M1 models inputs are already built here.

## Inputs

The relaxed VASP structures are kept as references in `CONTCARs`:

- `Cs2_Sn_Br6.CONTCAR`
- `Cs2_Sn_I6.CONTCAR`
- `Cs2_Ti_Br6.CONTCAR`
- `Cs2_Ti_I6.CONTCAR`
- `Cs2_Zr_Br6.CONTCAR`
- `Cs2_Zr_I6.CONTCAR`

The PySCF model inputs are in `Models`.

M0 folders contain:

- `quantum_region.xyz`

M1 folders contain:

- `quantum_region.xyz`
- `point_charges.csv`

## Models

- M0: isolated `[BX6]2-` octahedron with 7 QM atoms and no point charges.
- M1: `[Cs8BX6]6+` cluster with Cs atoms embedded in point charges.

The materials are:

- `Cs2_Sn_Br6`
- `Cs2_Sn_I6`
- `Cs2_Ti_Br6`
- `Cs2_Ti_I6`
- `Cs2_Zr_Br6`
- `Cs2_Zr_I6`

## Ground-State SCF

Run one ground-state DF-RHF calculation:

    python run_pyscf.py Models/Cs2_Sn_Br6/M0

Each successful SCF calculation writes:

- `scf.log`
- `scf.chk`

The log should contain `Status: CONVERGED`.

## Excited States and NTOs

Run one excited state calculation after SCF finishes:

    python run_tda.py Models/Cs2_Sn_Br6/M0


Each successful TDA calculation writes:

- `tda.log`

`run_tda.py` computes five singlet states and writes selected NTO information for State 1 and the brightest state among the first five states.

## Summary and Figures

After all `tda.log` files exist, make the summary tables:

    python summarize_tda.py Models

This writes:

- `Models/tda_summary.csv`
- `Models/tda_summary.txt`

Make the report figures:

    python make_figures.py

This writes:

- `Figures/figure_1_M0_M1_state1_energy.png`
- `Figures/figure_2_M1_state1_vs_brightest.png`

## Methods Summary

The electronic-structure workflow is:

1. DF-RHF ground-state calculation with `def2-SVP`.
2. ECPs for Cs, Zr, Sn, and I.
3. Five singlet TDA-HF/CIS excited states.
4. Selected NTO analysis for State 1 and the brightest root.
5. Summary table and basic report figures.
