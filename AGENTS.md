> **Before starting work in this directory, read [`../AGENTS.md`](../AGENTS.md)**
> for multi-agent coordination rules, build/test infrastructure, GHA workflows,
> and worktree discipline. That file is the authoritative reference for all
> cross-package agent operations.

# split-support: Agent Memory

## Project Overview

R-based simulation study evaluating character concordance measures and split
support metrics in phylogenetic inference. Accompanies the paper "Which
characters support which clades? Exploring the distribution of phylogenetic
signal using concordant information" (Smith, forthcoming).

## Directory Structure

```
split-support/
├── data-raw/
│   ├── _config.R                 # Central config: parameters, paths, helper fns
│   ├── 10_simulate.R             # Simulate alignments from reference tree
│   ├── 20_MrBayes*.R             # Bayesian inference (local + HPC variants)
│   ├── 30_iqtree.R               # Maximum-likelihood inference (IQ-TREE)
│   ├── 40_tnt.R                  # Parsimony inference (TNT)
│   ├── 80_byEdge.R               # Edge-wise analysis → Fig 4, Fig 3
│   ├── 90_byChar.R               # Character-wise analysis → Fig 5
│   ├── reference-gam.tre         # True reference tree (48 tips)
│   ├── alignments/               # 1 000 simulated NEXUS alignments
│   ├── MrBayes/                  # Bayesian results
│   ├── iqtree/                   # IQ-TREE results
│   ├── tnt/                      # TNT results
│   ├── concordance/              # Cached ClusteringConcordance results
│   ├── entropy/                  # Cached ClusteringEntropy results
│   ├── roc/                      # Cached per-panel roc + cIdx for Fig 4
│   └── nid/                      # Cached per-metric C-indices + GAM curves
│                                 #   for Fig 3 (see `nidCacheVersion`)
└── AGENTS.md
```

## Simulation Parameters (_config.R)

| Parameter | Value |
|-----------|-------|
| `nAln`    | 1 000 replicates |
| `nTip`    | 48 taxa |
| `nChar`   | 96 characters |
| `nCats`   | 6 gamma rate categories |
| `sim`     | `"gam"` (identifier used in all file-naming helpers) |

### File-naming helpers (from `_config.R`)

- `DataFile(sim, aln)` → NEXUS alignment path
- `MBFile(sim, aln, ext)` → MrBayes output path
- `IQFile(sim, aln, ext)` → IQ-TREE output path
- `TNTFile(sim, aln, type)` → TNT output path
- `ConcFile(sim, aln)` → concordance cache path
- `EntropyFile(sim, aln)` → entropy cache path
- `PartQFile(sim, aln)` → partition quality cache path

## Session Variables (populated by 80_byEdge.R)

| Variable | Type | Description |
|----------|------|-------------|
| `partCorrect` | logical vector | Whether each partition is in the reference tree |
| `partQual` | numeric vector | Normalised MCI of each partition vs. reference splits; 1 = true split |
| `postProb` | numeric vector | MrBayes posterior probabilities |
| `concord` | matrix | Concordance metrics: `quartet`, `wQuartet`, `cluster`, `phylo`, `mutual`, `shared`, `clusterNorm`, `wNrqs` |
| `tntStat` | matrix | TNT support: `symFq`, `symGC`, `boot`, `jak`, `pois` |
| `iqStat` | matrix | IQ-TREE support: `alrt`, `lbp`, `abayes`, `ufb` |
| `bremer` | numeric vector | Bremer (decay) support from TNT |
| `splitH` | matrix | Clustering entropy per partition |
| `allDat` | data.frame | All of the above combined (354 638 rows × 19 cols); includes NAs for method-specific metrics |
| `dat` | data.frame | `allDat` filtered to rows with no NAs in concordance + postProb columns |
| `refSplits` | Splits | Splits of the reference tree |
| `referenceTree` | phylo | The true 48-tip reference tree |

### Key derived quantity

```r
nid <- 1 - partQual  # Normalised clustering information distance (NID)
                     # 0 for every true split; (0, ~0.88] for incorrect splits
```

## Key Analyses (80_byEdge.R)

### Fig 4 — mosaic (spineplot) panels
Function `Histy(var, breaks, even, cf)` produces one mosaic panel showing how
`var` predicts `partCorrect`. Metrics are grouped by inference method (a =
MrBayes, b = IQ-TREE, c = TNT); concordance metrics appear in every group.
ROC-AUC and C-index are annotated on each panel. Output:
`Fig 4 - edge concordance.pdf` (5.4 × 8.4 in).

**Support chips.** The C-index and AUC are each printed on a rounded chip whose
fill encodes the value (`SupportAnnotation()` → `ChipRow()` → `RoundRect()`, all
inside the `# --- support chips:` block). Conventions:

- Ramp is **cividis** (`chipRamp`), chosen because it is monotonic in lightness
  (brighter = better reads directly) and its blue→yellow path avoids the
  red/green axis already spent on `partCorrect`. Ink switches between black and
  white by WCAG relative luminance (`ChipInk()`, threshold 0.179).
- The two statistics are scaled **independently**: `cIndexRange = c(0.5, 0.9)`,
  `aucRange = c(0.5, 1)`. Both are anchored at 0.5 (no discrimination) and
  stretched to their own attainable ceiling — the C-index never exceeds 0.86
  here. **Equal fills therefore do not mean equal values across the two rows**;
  `ChipLegend()` draws one bar per statistic in layout slot 25 (bottom right),
  in panel order (AUC, then C-index), under a short "better →" arrow.
- Geometry is computed in **inches** (`grconvertX/Y`) so corner radii stay
  isotropic under each panel's user coordinates. `ChipCex()` divides out the
  `par("cex") = 0.66` that `layout()` imposes, because `text()` scales `cex` by
  it whereas the `mtext()` call this replaced did not — without that the
  annotation would render at 4 pt rather than 6 pt.
- Panel titles are drawn by `SupportAnnotation()` at `line = 2.6`, *not* by
  `spineplot(main =)`, whose default placement leaves no room for the chips.
- **Both rows carry a 95% interval** at 3 dp, via `Bounds()`: `AUC = 0.79
  (0.784–0.789)` over `C-idx = 0.72 (0.723–0.725)`. AUC is the upper row, to
  match the order the manuscript introduces the two statistics in; `ChipLegend()`
  stacks its bars the same way. Two decimals would print bounds equal to the
  estimate — every interval here is narrower than the
  second decimal place (see the C-index section below). The C-index prefix is
  abbreviated to **`C-idx`** to buy the width: the widest row is then 1.05 in
  against the 1.35 in each panel's figure region allows, measured with
  `scratchpad/chip_mock.R`. `ChipLegend()` still spells out "C-index".
  Fig 3 deliberately shows **no** intervals — its panels have no room.

**Fig 4 metric order** (unique, first appearance — also used for Fig 3):
`cluster`, `mutual`, `wNrqs`, `wQuartet`, `postProb`, `ufb`, `lbp`, `abayes`,
`alrt`, `sCF`, `jak`, `boot`, `symFq`, `symGC`, `pois`, `bremer`

Exactly **two** quartet measures are plotted, both weighted, so that the
contrast is single-axis: `wNrqs` (non-redundant quartet statements,
chance-corrected — the new package default) comes **first**, ahead of
`wQuartet` (naive quartet currency, uncorrected — the old package default).
The unweighted `quartet` column is still computed and cached, but no longer
plotted. IQ-TREE's own `sCF` is retained as an independent measure, not as a
third variant of these.

### Fig 3 — NID vs. support scatter (referee request)
Function `.NidPanel(values, name)` plots NID (x) against a support metric (y).
Points are coloured by which methods carry data for the edge (`hereFor`); the
two trend lines are GAM fits, over `common` edges and over all edges with data.
Three C-indices are chipped in each panel — no intervals, for want of space.
Layout: 6 rows × 3 cols, 17 metric panels + 1 legend slot. Output:
`Fig 3 - CID vs support.pdf` (7 × 9 in). The stale `Fig A1 - NID vs
support.pdf` in the repo root predates the renumbering; nothing writes it.

**Support chips** (`.NidChips()`, same machinery as Fig 4 — `ChipRow()`,
`ChipFill()`, `ChipInk()`, `RoundRect()`). Three chipped rows down the
top-right corner *inside* each panel, replacing the plain `text()` annotations;
`nidChipCex = 0.78`, matching what that text used. Each of the three C-indices
gets its own bounds on the shared cividis ramp, all anchored on 0.5:

| Row | Subset | Domain | Observed |
|-----|--------|--------|----------|
| `C-index` | all edges with data for the measure | `cIndexRange` = 0.5–0.9 | 0.62–0.86 |
| `C_CTBI` | `common` edges only | `nidCtbiRange` = 0.5–0.95 | 0.56–0.93 |
| `C_CTBI'` | `common & !partCorrect` | `nidCtbi1Range` = 0.3–0.7 | 0.38–0.67 |

- The C-index row shares Fig 4's domain **deliberately**, so a chip means the
  same thing in both figures.
- `C_CTBI'` is the only one **symmetric** about 0.5, because it is the only one
  that goes below it (sCF = 0.38 — it ranks incorrect edges in the *wrong*
  order). Symmetry buys an absolute reading: mid-grey = uninformative, brighter
  = informative, darker = inverted. This is what makes the figure's point
  legible at a glance — concordance measures are olive/yellow on that row
  (0.56–0.67) while every resampling and likelihood measure is flat grey
  (0.47–0.51).
- Bounds are **hard-coded**, not derived at draw time, so the key reads the same
  between runs; `.NidChips()` warns on any value outside its domain, which
  `ChipFill()` would otherwise clamp silently.
- `ChipBar(labelSide = "left")` draws the three key bars into the legend slot
  beneath the "Stats available" point key, whose `cex` dropped 1 → 0.9 to make
  room.
- The chips are opaque and sit over the data, as the text they replaced did.

Values for the pre-chip version of this figure can be recovered from the text
layer of the committed PDF (cairo writes literal strings), which is how the
domains above were fixed without re-running the pipeline.

## Concordance Metrics (TreeSearch package)

| Column | Function | Notes |
|--------|----------|-------|
| `cluster` | `ClusteringConcordance(chanceCorrect = FALSE)` | Primary metric |
| `mutual` | `MutualClusteringConcordance()` | Primary metric |
| `wQuartet` | `QuartetConcordance(unit = "quartet", chanceCorrect = FALSE)` | Naive quartet currency |
| `wNrqs` | `QuartetConcordance(unit = "nrqs", chanceCorrect = TRUE)` | May be negative |
| `quartet` | as `wQuartet`, `weight = FALSE` | Cached; no longer plotted |
| `clusterNorm` | `ClusteringConcordance(chanceCorrect = TRUE)` | Not used in study |
| `phylo` | `PhylogeneticConcordance()` | Not used in study |
| `shared` | `SharedPhylogeneticConcordance()` | Not used in study |

## Key Dependencies

**R packages:** `TreeTools`, `TreeSearch`, `TreeDist`, `phangorn`, `ape`,
`survival` (C-index via `concordance`), `mgcv`, `pROC`, `viridisLite` (support
chips). **`Hmisc` is no longer used** — see the C-index note below.

**External software:** TNT, MrBayes 3.2.7, IQ-TREE 3.0.1

## The C-index (`SomersD` / `CIndex` in 80_byEdge.R)

Computed with **`survival::concordance()`**, not `Hmisc::rcorr.cens()`.

- **Why.** `rcorr.cens` is O(n²) (measured exponent 2.00): 19 min at
  n = 354 638. `.NidPanel()` needs 51 C-indices, six of them at full size, so
  Fig 3 alone took ~2 h per run. `survival::concordance()` uses an
  O(n log n) balanced tree — 2 s at the same n — and returns the identical
  statistic (verified to ≤1e-8 on continuous, integer-valued, heavily tied and
  negative-valued predictors against a `partQual`-shaped target).
- **Standard errors were wrong until 2026-08-05.** `rcorr.cens`'s `"S.D."` is
  already the standard *error* of Dxy; the old code divided it by `sqrt(n)`
  again, understating every interval by a factor of √n (595× at full n). The
  `ci95` element was never plotted, so no figure was affected, but **any C-index
  CI quoted in the manuscript from before this date is wrong.** Correct values:

  | Panel group | *n* | SE | 95% half-width |
  |---|---|---|---|
  | TNT (`tntCF`) | 36 761 | 0.0023–0.0044 | ±0.0045–0.0085 |
  | IQ-TREE (`iqCF`) | 44 955 | 0.0017–0.0033 | ±0.0033–0.0065 |
  | All-edge (`mlCF`) | 354 322 | 0.00047–0.00060 | ±0.0009–0.0012 |

  Widest of all 27: `wQuartet` in the TNT group, C = 0.568 ± 0.0085. Every
  interval is narrower than the second decimal place, which is why **Fig 4
  prints its bounds at 3 dp** (2 dp would repeat the estimate) and why Fig 3,
  which has no room for them, omits them entirely.
- `CIndex()` returns `estimate`, `ci95` **and `se`**. Do not `lapply` over it:
  `se` takes the scale change from Dxy to C but not the location shift.
- The `roc/` cache guard treats a stored `cIdx` with no `se` element as stale
  and rewrites it, so pre-fix caches heal themselves on the next run.

## Notes

- Concordance results are cached in `concordance/` and `entropy/` to avoid
  recomputation; delete cache files if partitions change.
- `data-raw/nid/` caches each Fig 3 metric's three C-indices *and* its two GAM
  prediction curves. The guard checks `nidCacheVersion` and nothing else, so
  **bump that constant** after any change to the GAM family, `gamma`, the
  spline basis, the subsets the C-indices use, or the meaning of a metric
  column. A second run of the figure then costs seconds, not hours.
- **Never call `QuartetConcordance()` or `ClusteringConcordance()` without
  pinning `unit` and `chanceCorrect`.** TreeSearch ≥ 2.0.0 renamed `normalize`
  to `chanceCorrect`, added `unit`, and defaults to
  `unit = "nrqs", chanceCorrect = TRUE`. A bare call therefore returns a
  *different measure* from the one the caches hold, and the cache guards in
  `80_byEdge.R` check only column names, never their meaning. Use the
  `NaiveQC()` / `NrqsQC()` helpers defined at the top of that script.
- `partQual = 1` for **all** splits in the reference tree by construction
  (they are assigned 1 before the MCI calculation loop).
- The `common` logical vector flags rows where all methods have non-NA values
  (~37–45k rows); use this for fair cross-method comparisons.
- Systematic Biology figure guidelines require `frame.plot = FALSE` in base R
  plots (no surrounding box).
