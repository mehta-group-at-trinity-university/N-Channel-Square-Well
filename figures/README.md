# Figure and analysis notebooks

Mathematica notebooks that read the Fortran output in `../paper-data/` and produce the figures and
statistics for Phys. Rev. A **98**, 062703 (2018). Recovered 2026-09-22 from
`~/Documents/Code/Nchannel/`, where they had never been under version control.

## Running them: the path problem

The `Import` statements use **absolute paths** from an older layout, e.g.

    Import["~/Code/Nchannel/Square/CODES/resdata_from_morgoth/fort.800", "Table"]

That tree no longer exists (it became `~/Documents/Code/`, now being retired). The data those
statements want is in `../paper-data/`, laid out with the **same relative structure**, so the
simplest way to run a notebook unchanged is a symlink:

```bash
mkdir -p ~/Code
ln -s ~/Developer/N-Channel-Square/paper-data ~/Code/Nchannel
```

Then `~/Code/Nchannel/Square/CODES/...` resolves into this repo. Alternatively, edit the `Import`
paths (note that the `.nb` format wraps strings at 80 columns, so a filename may be split across
lines — edit in Mathematica, not a text editor).

## Which notebook needs which data

| Notebook | Imports |
|---|---|
| `makeScatStatPlots*.nb` (4 versions) | `Square/CODES/fort.{500,700,800}`, `resdata_from_morgoth/fort.{500,700,800}`, `resdata_from_morgoth/end-cases/fort.800` |
| `makeScatStatPlotsv3/v4` additionally | `Square/CODES/plots/linespacings.dat`, `plots/NETable.dat` |
| `MakeChrisTicnorPlots*.nb` | `Square/CODES/fort.800`, `fort.100` |
| `widthstats.nb` | `Square/CODES/fort.800` |
| `BrodyStatisticsV2.nb` (not in this repo, see below) | `fort.{800,801,700,1001}`, `gaussian_voffdiag_from_morgoth/fort.{500,700}` |

`fort.801` and `MagneticF/fort.100` were recovered from inside `Square/NChannelCodes.tar`; the
others were loose files.

## Two gaps

1. **`gaussian_voffdiag_from_morgoth/`** — imported by `BrodyStatisticsV2.nb`, but **this data no
   longer exists anywhere**: not in the old tree, not in any tarball, not in this repo. It was
   presumably left on the cluster (morgoth). Those cells cannot be re-run. The *results* survive
   only as stored output inside the notebook, which is the main reason that file is 182 MB — so do
   not clear its output.
2. **`BrodyStatisticsV2.nb` (182 MB) is not in this repo.** It exceeds GitHub's 100 MB per-file
   limit. It is preserved in `~/Documents/Code-source-archive-2026-09-22.zip`. To version it,
   either track it with Git LFS or store it outside git. The smaller `BrodyStatistics.nb` (42 MB,
   2018-05-28) *is* here, in `../bound-sector/`.

## Naming

Files ending `-legacy.nb` are the `~/Documents/Code` copies of notebooks whose same-named versions
already existed at the repo root and **differ from them** (`MakeChrisTicnorPlots`,
`makeScatStatPlotsv3`, `Time-Delay-Formulas`). Both were kept rather than guessing which is
authoritative — worth reconciling when you next open them.

## Versions

`makeScatStatPlots` exists as four progressive versions (base 2018-08-16, LogScale 2018-07-25,
v2 2018-08-28, v3 2018-09-14, v4-no-widths 2018-11-16), and `CoupledSqrWell` / `MyNchannelSqrWell`
similarly. All are kept; the highest version number is usually, but not reliably, the one used for
the paper — the published figures are dated 2018-08-17.
