# Bound-state sector

Multichannel **bound-state** solver for the same model whose scattering sector is computed by
`../squarenew.f90` (published as Phys. Rev. A **98**, 062703 (2018)).

From the header of `MCBound.f90`:

> THIS PROGRAM CALCULATES THE BOUND STATES SUPPORTED BY A MULTICHANNEL SET OF POTENTIALS THAT ARE
> SPECIFIED BY A SUBROUTINE. NB: Because this is a bound-state code, it does NOT do R-matrix
> propagation.

Method: B-spline basis, overlap and Hamiltonian matrices assembled by Gauss-Legendre quadrature,
solved as a generalized eigenvalue problem through `Mydggev` (LAPACK `DGGEV`).

## Why this belongs with `squarenew.f90`

`MCBound.f90`'s `makeVTridiag` is **byte-identical** to the one in `../squarenew.f90` — the same
tridiagonal coupling model, in closed form. This is the bound-state counterpart of the published
scattering calculation, not a separate model.

## Relationship to `LoopOverCalcBound.f`

Two variants of the same solver exist, differing only in where the potentials come from:

| | this directory | `~/Developer/LiScattering/LoopOverCalcBound.f` |
|---|---|---|
| Potentials | **closed form** via `makeVTridiag` | **interpolated adiabatic curves** via `ReadCurves`, `setup_interp`, `SetMultipoleCoup` |
| Pipeline | matches the square-well paper | CABS-style (adiabatic U/P/Q) |
| Newest version | `MCBound.f90`, 2020-07-03 | 2021-06-19 |

They share 12 routines: `CalcBoundStates`, `CalcHamiltonian`, `CalcOverlap`, `CheckBasis`,
`CheckBasisPhi`, `CheckPot`, `CompSqrMatInv`, `deigsrt`, `GridMaker`, `Mydggev`, `printmatrix`,
`setup_potential_matrix`. Neither supersedes the other — if you change one, consider the other.

The 2013 `LoopOverCalcBound.f` kept here (from `MCBoundStates.working.tar`) is the common
ancestor, 93 lines different from LiScattering's maintained copy. A third copy lives in
`~/Developer/AtomIon1D/MCBoundStates.working/`.

## Contents

| File | Date | Notes |
|---|---|---|
| `MCBound.f90` | 2020-07-03 | The solver. Newest source in the original `Nchannel` tree |
| `MCBound.par`, `MCBound.mak` | 2018-06-16 | Parameters and build file |
| `LoopOverCalcBound.f`, `.mak`, `calcbound.par` | 2013-06 | Ancestor, interpolated-curve line |
| `BoundStateCalc.nb` | 2017-07-14 | Bound-state analysis |
| `MySparseBound.nb` | 2017-06-27 | Sparse-matrix approach |
| `BrodyStatistics.nb` | 2018-05-28 | Level-statistics analysis (42 MB) |

Dependencies (`Bsplines.f`, `matrix_stuff.f`, `besselnew.f`, `modules_qd.f90`) are not duplicated
here — they are byte-identical to copies in `~/Developer/lib` and elsewhere in this repo.

## Provenance

Recovered 2026-09-22 from `~/Documents/Code/Nchannel/BoundSector/`, which was otherwise build
artifacts (`.o`, `.mod`, `.x`) and `fort.*` output. Never previously under version control.

## Next steps

The multi-open-channel generalization of this model is being developed in
[`Erbium`](https://github.com/quantumwave/Erbium) (Mathematica). The Fortran here, like
`squarenew.f90`, currently assumes a single open channel.
