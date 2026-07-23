# FunTable

**FunTable** is a Wolfram Language (Mathematica) package for generating numerical tables of
**conformal-bootstrap functionals** (1d / 2d / 3d) and assembling them into linear-programming
(LP) bootstrap problems.

It evaluates the 1d crossing functionals in arbitrary precision, builds the 2d and 3d
functionals from them by dimensional reduction, and exports packed **HDF5** tables for
downstream LP solvers.

## Requirements

- Wolfram Language / Mathematica (developed with 13.x–14.x).
- A working C compiler (`gcc`/`clang`). Performance-critical kernels in `src/Table.wl` are
  compiled with `CompilationTarget -> "C"` and cached in `lib/` on first load; the cache is
  platform-specific and is (re)built automatically if absent.

## Loading

Load the package **as a file**, so `$InputFileName` resolves the `src/` paths:

```mathematica
Get["/path/to/FunTable/FunTable.wl"]
```

`FunTable.wl` sets its own directory and loads `src/OneDfunctional.wl`, `src/OneDAsymp.wl`,
and `src/Table.wl` in order.

## Quick start

Generate a 3d functional table (HDF5) at external dimension Δφ:

```mathematica
Get["/path/to/FunTable/FunTable.wl"];   (* as a file *)
LaunchKernels[];                        (* the 1d generation step is parallel *)

dphi3d = 5181488023/10000000000;        (* exact rational; e.g. the 3d Ising point *)
dphi1d = dphi3d/2;                      (* convention: 1d dphi = 3d Δφ / 2 *)

(* 1. 1d table (.mx): orderMax = 10, ΔMax_gen = 650, exp-normalized *)
generate1ptmx["Ising.mx", dphi1d, 10, 650, True];

(* 2. 3d table: normalized, with 30 sub-bound "ghost" interpolation rows *)
lp = init3DTablePatch["Ising.mx", 10, Range[0, 100, 2], 150, 500,
                      "normed" -> True, "fill" -> 30];

(* 3. export to HDF5 *)
lpN = KeyDrop[KeyMap[# /. "Δtable" -> "dtable" &, lp], "normed"];
Export["Ising.h5", Normal[lpN], "HDF5"];
```

The **`FunTable.nb`** notebook is a runnable version of this (1d single-point evaluation +
the full 3d-table pipeline). For the complete technical reference see
**[`H5_generation_technical_summary.md`](H5_generation_technical_summary.md)**.

## Repository layout

| path | contents |
|---|---|
| `FunTable.wl` | top-level loader; defines `generate1ptmx` |
| `src/OneDfunctional.wl` | 1d crossing functionals (arbitrary-precision, hypergeometric) |
| `src/OneDAsymp.wl` | large-Δ asymptotic expansions of the same functionals |
| `src/asympcoeff.wl` | hardcoded asymptotic-series coefficients |
| `src/Table.wl` | 2d/3d construction, normalization, HDF5 tables, LP setup |
| `FunTable.nb` | usage-example notebook |
| `H5_generation_technical_summary.md` | technical reference for the generation pipeline |
| `doc/` | focused notes (exp-only normalization, `init3DTablePatch`) |
| `lib/` | compiled-kernel cache (generated on first load; not shipped) |

## Key entry points

- `generate1ptmx[file, dphi1d, orderMax, ΔMax_gen, normed]` — build and save the 1d table (`.mx`).
- `init1DTable` / `init2DTable` / `init3DTable` / `init3DTablePatch` — load a saved table or
  derive a higher-dimensional one; return an `Association` (`"functablePacked"`, `"dtable"`,
  `"spinlist"`, `"order"`, …).
- `sectionPosList`, `gapcostlist`, `initLP` — build the Δ index slabs and the LP problem.

## Conventions

- **1d dphi = 3d Δφ / 2** (external dimension).
- Δ grids: 1d `Range[0, ΔMax_gen, 1/202]`; 2d/3d spacing `1/101` (`halfstep = 101`), index
  `i` ↔ `Δ = (i−1)/101`.
- Use **exact rationals** for Δφ so the high-precision generation is not limited by float
  rounding of the input.
- Generation configs: **High** (`ΔMax_gen = 500`, `nmax = 350`), **VeryHigh** (`650`, `500`);
  `orderMax` is the functional truncation order (10 standard, 12 extended).

Normalization, the ghost-node `fill`, the packed data layout, and special-Δφ handling are
documented in [`H5_generation_technical_summary.md`](H5_generation_technical_summary.md).

## License

*(No license file yet — add one before public release.)*
