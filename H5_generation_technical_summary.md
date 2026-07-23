# Generation of 3d conformal-bootstrap functional tables (`.h5`)

Technical summary of the `FunTable` pipeline that produces the packed 3d crossing-functional
tables (HDF5) consumed by the linear-programming (LP) bootstrap solver. Symbol names and
file references point at the source (`src/OneDfunctional.wl`, `src/OneDAsymp.wl`, `src/Table.wl`,
`FunTable.wl`).

---

## 0. Overview

The 3d functionals are built by **dimensional reduction from 2d**, and the 2d functionals from
**1d crossing functionals**. The pipeline is two stages:

1. **1d table** (`generate1ptmx` → `.mx`): evaluate the 1d crossing functionals on a uniform
   Δ grid, stitching a *direct* (arbitrary-precision hypergeometric) low-Δ region to an
   *asymptotic* (large-Δ series) high-Δ region, and `DumpSave` the result.
2. **3d table** (`init3DTablePatch` → `.h5`): from the 1d table build the 2d functionals
   (`func2dspin`), then the 3d functionals per spin (`func3dspin`) as a truncated descendant
   sum over 2d blocks plus an analytic tail, apply the standard normalization, and `Export`
   to HDF5.

**Conventions.** The 1d external dimension is half the 3d one, `dphi₁ᴅ = Δφ₃ᴅ / 2`. The 1d Δ
grid is `Range[0, ΔMax_gen, 1/202]` (spacing `1/202`); the 2d/3d grid is `Δtable2 = 2·Δtable1`,
spacing `1/101`. The integer `halfstep = 101` encodes this spacing: grid index `i` ↔
`Δ = (i−1)/101`. Use **exact rationals** for `dphi` so the high-precision generation is not
limited by float rounding of the input.

---

## 1. 1d crossing functionals (`src/OneDfunctional.wl`)

Four functionals are tabulated: the **bosonic ("minus")** pair `betaminus`, `alphaminus` and
the **fermionic ("plus")** pair `betaplus`, `alphaplus`. Each is evaluated at external
dimension `Δφ`, functional order `m`, and operator dimension `Δ`.

**Three-term recursion.** The building blocks `b0[n]`, `a0[n]` (minus) and `b1[n]`, `a1[n]`
(plus) satisfy a linear three-term recursion with coefficients

- `μ[n] = −n²`,  `μ′[n] = −2n`
- `ν[Δφ,Δ][n] = (−1+Δ)Δ + Δφ + ½ n(−1+n+4Δφ)`
- `ρ[Δφ,Δ][n] = −((n+2Δφ)² (−1+n+4Δφ)²) / (4(−1+2n+4Δφ)(1+2n+4Δφ))`
- `ν′`, `ρ′` (the α–β coupling terms), with
  `ρ′[Δφ,Δ][n] = −((n+2Δφ)(−1+n+4Δφ)(1+4n³−6Δφ+24n²Δφ+32Δφ³+n(−2+48Δφ²))) / (2(−1+2n+4Δφ)²(1+2n+4Δφ)²)`

seeded by `b0[−1]=0`, `b0[0]=β0s0OverA` (a closed-form `HypergeometricPFQRegularized`), and
driven by inhomogeneous terms `ROverA`, `SOverA`:

```
b0[n] = ( ROverA[(n−1)/2] − ( ρ[n−2]·b0[n−2] + ν[n−1]·b0[n−1] ) ) / μ[n]
a0[n] = ( SOverA[(n−1)/2] − ( ρ[n−2]·a0[n−2] + ν[n−1]·a0[n−1]
                              + ρ′[n−2]·b0[n−2] + ν′[n−1]·b0[n−1] + μ′[n]·b0[n] ) ) / μ[n]
```

The functionals are the "self-shadow-symmetric" aggregates
`βSS[Δφ,Δ][n] = b0[2n] + 2·bCS[n] + bs0OverA[n]` (and `αSS` analogously), returned negated:
`betaminus[Δφ,m,Δ] = −βSS[…][m]`, etc. `bCS`, `bs0OverA` are closed-form correction terms.

**Aggregators.**
- `bosonicminusfuncList[Δφ,order,Δ]` → `{ {0}~Join~(−βSS/@Range[order]), −αSS/@Range[0,order] }`
  — a `2 × (order+1)` list (β-row padded with a leading 0, α-row full).
- `fermionicplusfuncList[Δφ,order,Δ]` → `{ −βpF/@Range[0,order], −αpF/@Range[0,order] }`.
- `functionalpair[Δφ,order,Δ] = { bosonicminusfuncList, fermionicplusfuncList }`
  (shape `2 × 2 × (order+1)`).

**Precision.** Evaluation is arbitrary-precision (`HypergeometricPFQRegularized`,
`SetPrecision` to the `PrecisionGoal` option, `defaultPrec = 150`; the α-functionals work at
`3/2 × PrecisionGoal`). Both `OneDfunctional`` and `OneDAsymp`` register in
`$DistributedContexts` so the symbols are available on parallel kernels.

**Degenerate-point fix.** `ρ` and `ρ′` carry a factor `(−3+4Δφ)` in the `n=−1` denominator,
which vanishes at `Δφ = 3/4`, turning the (formally zero) `ρ[−1]·b[−1] = ∞·0` into
`Indeterminate` and poisoning the whole direct region. Since these `n=−1` terms always multiply
the seed `b[−1]=0`/`a[−1]=0`, they must be exactly zero; the fix defines
`ρ[Δφ_,Δ_][−1] = 0` and `ρ′[Δφ_,Δ_][−1] = 0` (identical fix mirrored in `OneDAsymp.wl`).

---

## 2. Asymptotic (large-Δ) branch (`src/OneDAsymp.wl`)

Direct hypergeometric evaluation is prohibitively slow at large Δ, so the high-Δ region uses
**asymptotic expansions** of the same functionals: `bosonicminusAsymp`, `fermionicplusAsymp`
(`qmax = 27+2·order`, working precision `prec = 2·qmax−10`). They use the same recursion with
Δφ-only coefficients and carry the same `ρ[−1]=ρ′[−1]=0` fix.

- `asymp1dfast` / `asymp2dfast` return the 2d asymptotic **series pair** `(serminus, serplus)`
  as `SeriesData` truncated to `Δsermax` terms.
- `asymptablelist[Δφ,order,Δlist]` produces the `{minushigh, plushigh}` blocks that fill the
  high-Δ part of the 1d table.
- `series2value[ser,Δ] = Σᵢ ser[[i]] · Δ^(−5−i)` evaluates such a series at a given Δ.
- `functionallistCutOff[Λ]` enumerates the independent functional components at truncation
  order `Λ` (built from the base pairs `{{1,1},{2,1},{1,2},{2,2}}` and the cut
  `functionalindCut`); its length is the exported `funcnum` (e.g. 242 at order 10, 338 at
  order 12).
- The large-Δ resummation uses **hardcoded extended series** (the `AnjlD3`/`Aser` block and
  the `nplus`/`nminus` normalization ratios), computed out to `kmax = 16`
  (`data/aser_kmax16/`, `data/normfac_kmax16/`).

---

## 3. 1d table generation — `generate1ptmx` (`FunTable.wl`)

```
generate1ptmx[filename, dphi₁, orderMax, ΔMax_gen, normedFlag]
```

- **Grid & split.** `table = Range[0, ΔMax_gen, 1/202]`. A threshold
  `lowseuil = 2·dphi + 2·orderMax + 25 + EulerGamma` separates `tablelow` (direct evaluation)
  from `tablehigh` (asymptotic).
- **Dynamic precision.** `prec` is raised in steps of 10 until a sample evaluation reaches
  ≥ 27 significant digits: `While[precMax < 27, prec += 10; precMax = Precision@functionalpair[dphi, orderMax, tablelow[[-1]], PrecisionGoal→prec]]`.
- **Parallel evaluation.** `minuslow`/`pluslow` via `ParallelMap[…, Method→"FinestGrained"]`
  over `tablelow` (the fermionic map is wrapped in `TimeConstrained[…,2]`); the high-Δ blocks
  come from `asymptablelist`; the two regions are `Join`ed into `minus`, `plus`.
- **Exp-only normalization** (when `normedFlag`): `minus = N[2^(−2·Δtable) · minus]` (same for
  `plus`), performed **in arbitrary precision before conversion to machine reals**. This is
  essential: the raw functionals grow like `4^Δ`, which overflows IEEE doubles (`~10³⁰⁸`) for
  `Δ ≳ 256`; rescaling by `2^(−2Δ)` while still at ~50-digit precision gives
  `2^(−2Δ)·4^Δ·poly = poly`, which `N[]` then converts safely.
- **Output.** `DumpSave[filename, {dphi, precMax, orderMax, Δtable, minus, plus, normed}]`.
  (`dphi` is stored as a machine real via `N@dphi₁`, even though the tables themselves were
  computed at the exact/high-precision input.)

---

## 4. 2d construction — `func2dspin` (`src/Table.wl`)

For spin `l`, the 2d functional is assembled from the shifted 1d data on the doubled grid
`Δtable2 = 2·Δtable1` (`halfstep = 101`):

```
res[[All, l·halfstep+1 ;; len − l·halfstep]] =
   ½ ( plus1[[All, 1 ;; len−2l·halfstep]] · minus1[[All, 2l·halfstep+1 ;; len]]
     + minus1[[All, 1 ;; len−2l·halfstep]] · plus1[[All, 2l·halfstep+1 ;; len]] )
```

i.e. a symmetrized product of the "plus" and "minus" 1d functionals shifted by ±`l` in Δ,
populated from index `l·halfstep+1` (`Δ = l`, the 2d unitarity edge) upward.

---

## 5. 3d construction — dimensional reduction (`func3dspin`, `src/Table.wl`)

The 3d functional of an operator of dimension Δ and spin `ℓ` is a weighted sum of 2d
functionals of the tower `(Δ+2n, 2j)`:

```
res[[funciter, Δiter]] = Σ_{jhalf=0}^{ℓ/2} Σ_{n=0}^{nmax}
     factor3dMatrix[jhalf, n] · listPacked2d[[funciter, jhalf+1, Δiter + n·2·halfstep]]
```

- **`factor3dMatrix`** is the cumulative product over descendants `n` of
  `factor3dCompiled[ℓ, n, jhalf, Δ]`, seeded (column `n=0`) by `ZListCompiled[ℓ]`.
  `factor3dCompiled` is a rational function of `(Δ, n, ℓ, jhalf)` with three branches
  (a `jhalf=0` case with an `n=1, ℓ=0` special sub-case; an `ℓ=2·jhalf` case; the general
  case). It carries a **twist-zero pole at `Δ = ℓ`** (a `1/(Δ−ℓ)` factor at `n=1`).
- **`ZListCompiled[ℓ]`** are the `j = 0,2,…,ℓ` prefactors,
  `(2/π)·Exp[LogGamma(½(−j+ℓ+1)) + LogGamma(½(j+ℓ+1)) − LogGamma(½(−j+ℓ+2)) − LogGamma(½(j+ℓ+2))]`,
  with the `j=0` entry halved.
- The compiled kernels (`func3dspin`, `func3dspinNormed`, `func3dspinAccumulate`,
  `func2dspin`, `factor3d*`, `ZListCompiled`, …) are cached in `lib/` as
  `<sym>.so` + `<sym>lib.mx` via `creatloadCompile` (built with `CompilationTarget→"C"`,
  auto-rebuilt if missing or if loading fails). `func3dspinNormed` inserts a factor `16` per
  descendant step; `func3dspinAccumulate` returns running descendant totals (an extra
  `nmax+1` axis).

**Truncation.** The direct descendant sum runs to `n = nmax`; the `n > nmax` remainder is
supplied by the analytic tail (§7).

---

## 6. Unitarity-bound layout and the ghost-node fill

`func3dspin` initializes `res` to all zeros and populates it only from index
`halfstep·(ℓ+1)+1 − fill`, i.e. from `Δ = ℓ+1 − fill/101`. In addition, the `ℓ = 0` case fills
index 1 (`Δ = 0`) with `listPacked2d[[funciter,1,1]]` (the scalar/identity contribution). So a
row for spin ℓ has the structure:

| region | content |
|---|---|
| `Δ = 0` | populated for `ℓ = 0` only (identity); zero for `ℓ ≥ 2` |
| `0 < Δ < ℓ+1 − fill/101` | zero (padding) |
| `Δ ≥ ℓ+1 − fill/101` | populated functional values |

- `fill = 0` reproduces the classic layout (populated only for `Δ ≥ ℓ+1`, the unitarity bound).
- `fill = 30` (the default) additionally fills the **30 grid points immediately below** each
  spin's bound with the analytic continuation of the descendant sum (same formula, evaluated
  at lower Δ). These "ghost" rows let near-bound interpolation use a **centered** stencil
  instead of a one-sided edge window; they are interpolation nodes only and never enter the LP
  as operators. The fill is safe up to `fill ≈ 71` before the lowest descendant read reaches
  the 2d table's `Δ = ℓ` edge (which is also the `factor3d` twist-zero pole); at `fill = 30`
  the lowest node sits `~0.70` in Δ above that pole.

---

## 7. The analytic tail beyond `nmax` (`asympspinlist`, `src/Table.wl`)

`asympspinlist` supplies the `n > nmax` remainder of the descendant sum via
`HurwitzZeta[−s, nmax+1+Δ/2]` (closed-form resummation of the descendant index) times the
asymptotic series (`Aser`/`normfac`). It is defined on `spanasympl[ℓ] = ((ℓ+1)·halfstep+1 ;;)`,
i.e. only for `Δ ≥ ℓ+1`, and is `PadLeft`-zeroed below the bound. Near the bound its magnitude
is tiny (`~10⁻¹⁹` for the tables here), so it is a small correction to the truncated direct sum.

---

## 8. Normalization

The final table is optionally multiplied by a spin- and Δ-dependent normalization
(`init3DTablePatch` option `"normed"→True`). Because the 1d `.mx` is stored already
exp-normalized (§3), the polynomial form is used:

```
normtoprodCompile3dPoly:  (1 + ½(Δ+ℓ))^{2Δφ+3/2} · (1 + ½(Δ−ℓ))^{2Δφ+7/2} / √(1+Δ)
```

populated from `Δ = ℓ` upward. (The un-normalized-mx path `normtoprodCompile3d` additionally
carries `4^{−Δ}`.) Here `Δφ` denotes the **3d** external dimension (`dphi = 2·lp1["dphi"]`).
The two factors `(1+½(Δ±ℓ))` correspond to the two shifted ("ρ"/"τ") 1d arguments, and the
`1/√(1+Δ)` is the extra `√Δ` factor carried by the 3d functional relative to the raw 1d one.

---

## 9. Assembly and HDF5 format — `init3DTablePatch` (`src/Table.wl`)

```
init3DTablePatch[lp1_or_file, order, spinlist, Δmax, nmax, "normed"→True, "fill"→30]
```

Steps: load/accept the 1d table `lp1`; select the functional components with
`functionallistCutOff[order]`; set `step = 202`, `halfstep = 101`, `Δtable2 = 2·Δtable1`; build
`listPacked2d` with `func2dspin` (chunked in groups of 30 spins); compute
`func3dspin`/`func3dspinNormed`, **truncate** the Δ axis to `Δmax·halfstep+1`; **add**
`asympspinlist`; apply the normalization. The returned Association (run scripts rename
`"Δtable"→"dtable"` and drop the `"normed"` key before `Export[…, "HDF5"]`) holds:

| key | meaning |
|---|---|
| `Dim` | `3` |
| `dphi` | 3d external dimension `Δφ` |
| `funcnum` | number of functional components = `Length@functionallistCutOff[order]` |
| `spinlist` | spins tabulated, e.g. `Range[0,100,2]` |
| `nmax` | descendant cutoff |
| `halfstep` | `101` (grid spacing `1/101`) |
| `dtable` | Δ grid `Range[0, Δmax, 1/101]`, length `Δmax·101 + 1` |
| `functablePacked` | the table, dims `{nspins, nΔ, funcnum}` |
| `order` | functional truncation order |
| `prec` | precision recorded in the 1d table |

`functablePacked[spin, Δ, component]`; the `Δ` axis follows the unitarity-bound layout of §6.

---

## 10. Generation configurations

| config | `orderMax`/`order` | `spinlist` | `ΔMax` (LP cutoff) | `ΔMax_gen` (1d) | `nmax` |
|---|---|---|---|---|---|
| **High** | 10 | `Range[0,100,2]` | 150 | 500 | 350 |
| **VeryHigh** | 10 | `Range[0,100,2]` | 150 | 650 | 500 |
| **VeryHigh, extended order** | 12 | `Range[0,100,2]` | 150 | 650 | 500 |

All with `normed→True` (and `fill→30` for the ghost-node tables). The 1d `.mx` generation is
the parallel step (`ParallelMap`); the 3d h5 conversion is largely single-kernel (compiled
`func3dspin`). Because the `.mx` is fill- and `nmax`-independent above the bound, a fixed 1d
table can be reused across `fill`/`order`/`nmax` variants of the 3d conversion.

---

## 11. Special-point handling

Certain 1d `dphi` require care and are worth recording for reproducibility:

1. **`dphi₁ᴅ = 3/4`  (`Δφ₃ᴅ = 3/2`)** — the `ρ`/`ρ′` `n=−1` denominator vanishes; fixed at
   source in both `OneDfunctional.wl` and `OneDAsymp.wl` (§1). Exact rationals then work.
2. **`dphi₁ᴅ = 1/2`  (`Δφ₃ᴅ = 1`)** — `fermionicplusfuncList` fails to return at small
   odd-integer Δ (degenerate hypergeometric parameters), in both exact and inexact arithmetic.
   Workaround: evaluate at `SetPrecision[dphi, 200] + 10⁻³⁰` — a `10⁻³⁰` shift that lifts the
   degeneracy; the induced table error `~10⁻³⁰` is far below the 27-digit target.
3. **General.** Values where `4·dphi` (or hypergeometric parameter combinations) hit
   (half-)integers can misbehave; check before a production run. Passing the external dimension
   as an **exact rational** (rather than a decimal) avoids float rounding limiting the
   high-precision evaluation.

---

## 12. Downstream LP construction (brief)

The table feeds the LP setup (`src/Table.wl`):

- `sectionPosList[Dim, Λ, dphi, halfstep, spinlist, Δtable2, Δcutoff, smallstep, bigstep]`
  builds the list of `{spin, Δ}` indices to `Extract` from `functablePacked`: a dense slab
  (stride `smallstep`) up to `sectionΛ[Λ,dphi,halfstep]` and a sparse slab (stride `bigstep`)
  out to `Δcutoff·halfstep+1`, each offset by `spin·halfstep`.
- `gapcostlist[lgap, LList, ΔList]` assigns the cost vector for a gap-maximization at spin
  `lgap`.
- `initLP[functionalTable, lgap, …]` extracts the constraint matrix `vecM` and target
  (`−functablePacked[[1,1]]`, the unit/identity row), sets up the auxiliary basis and effective
  costs (via `LinearSolve`), and returns the LP problem Association.

Sub-bound ghost rows (§6) participate only as interpolation nodes in this step; they are never
admitted as candidate operators.
