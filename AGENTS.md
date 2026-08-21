# split-support: Agent Memory

## Project Overview

R-based simulation study evaluating character concordance measures and split
support metrics in phylogenetic inference. 
Accompanies the paper "Which characters support which clades? Exploring the distribution of phylogenetic signal using mutual information" (Smith, in production, Systematic Biology).

## Directory Structure

```
split-support/
├── data-raw/
│   ├── _config.R                 # Central config: parameters, paths, helper fns
│   ├── 10_simulate.R             # Simulate alignments from reference tree
│   ├── 20_MrBayes*.R             # Bayesian inference (local + HPC variants)
│   ├── 30_iqtree.R               # Maximum-likelihood inference (IQ-TREE)
│   ├── 40_tnt.R                  # Parsimony inference (TNT)
│   ├── 80_byEdge.R               # Edge-wise analysis → Fig 3, Fig 4
│   ├── 90_byChar.R               # Character-wise analysis → Fig 5
│   ├── Figure_tree.R             # Reference topology → Fig 2
│   ├── reference-gam.tre         # True reference tree (48 tips)
│   ├── mb-gam.nex                # MrBayes block (brlenspr rate 0.245568)
│   ├── tnt-ew.run, bremer.run    # TNT scripts: search, support, Bremer
│   ├── slurm.sh                  # SLURM template for the HPC route
│   ├── alignments/               # 1 000 simulated NEXUS alignments
│   ├── MrBayes/                  # Bayesian results
│   ├── tnt/                      # TNT results
│   ├── iqtree/                   # IQ-TREE results (incl. `.cf.tree` sCF trees)
│   ├── concordance/              # Cached ClusteringConcordance results
│   ├── entropy/                  # Cached ClusteringEntropy results
│   ├── roc/                      # Cached per-panel roc + cIdx for Fig 4
│   └── nid/                      # Cached per-metric C-indices + GAM curves
│                                 #   for Fig 3 (see `nidCacheVersion`)
├── .zenodo.json                  # Metadata for the Zenodo deposit
├── README.md
└── AGENTS.md
```

`roc/` and `nid/` are **not under version control**: they hold derived results
that the figure scripts rebuild on demand, and `roc/` alone ran to 200 MB.
`_config.R` creates both, so a fresh clone works without them.
No figure is kept in the repository — each script writes its own, and all five
appear in the published paper. `trit-analysis/`, a standalone side-study that
feeds no published figure, is ignored too.

## Simulation Parameters (_config.R)

| Parameter | Value |
|-----------|-------|
| `nAln`    | 1 000 replicates |
| `nTip`    | 48 taxa |
| `nChar`   | 96 characters |
| `nCats`   | 6 gamma rate categories |
| `sim`     | `"gam"` (identifier used in all file-naming helpers) |

### File-naming helpers (from `_config.R`)

- `DataFile(sim, id, ext = ".nex")` → NEXUS alignment path
- `MBFile(sim, id = "", suffix = NULL)` → MrBayes output path
- `IQFile(sim, id, suffix = "")` → IQ-TREE output path
- `TNTFile(sim, id, wt = "ew")` → TNT output path
- `ConcFile(sim, id, suffix = "")` → concordance cache path; `90_byChar.R`
  passes `"_chr"`, `"_chrQ"` and `"_cns"` for its per-character caches
- `EntropyFile(sim, id)` → entropy cache path
- `PartQFile(sim, id)` → partition quality cache path

Directory constants: `alnDir`, `mbDir`, `iqDir`, `tntDir`, `concDir`, `hDir`,
`rocDir`, `nidDir`. **`tntDir` is `data-raw/tnt/`, lower case.** The case must
match the directory that ships with the repository: on a case-sensitive
filesystem a mismatch yields an empty directory, and every replicate is then
skipped for want of TNT results.

## Session Variables (populated by 80_byEdge.R)

| Variable | Type | Description |
|----------|------|-------------|
| `partCorrect` | logical vector | Whether each partition is in the reference tree |
| `partQual` | numeric vector | Normalised MCI of each partition vs. reference splits; 1 = true split |
| `postProb` | numeric vector | MrBayes posterior probabilities |
| `concord` | matrix | Concordance metrics: `quartet`, `wQuartet`, `cluster`, `phylo`, `mutual`, `shared`, `clusterNorm`, `wNrqs` |
| `tntStat` | matrix | TNT support: `symFq`, `symGC`, `boot`, `jak`, `pois` |
| `iqStat` | matrix | IQ-TREE support: `alrt`, `lbp`, `abayes`, `ufb`, `sCF` |
| `bremer` | numeric vector | Bremer (decay) support from TNT |
| `splitH` | matrix | Clustering entropy per partition |
| `allDat` | data.frame | All of the above combined (354 638 rows × 23 cols); includes NAs for method-specific metrics |
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
  `par("cex") = 0.66` that `layout()` imposes: `text()` multiplies its `cex` by
  it, so without that division the annotation renders at 4 pt rather than 6 pt.
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

**Fig 4 metric order** (unique, first appearance — 16 metrics over 24 panels,
since the four concordance measures repeat in each method group):
`cluster`, `mutual`, `wNrqs`, `wQuartet`, `postProb`, `ufb`, `lbp`, `abayes`,
`alrt`, `sCF`, `jak`, `boot`, `symFq`, `symGC`, `pois`, `bremer`

Fig 3 follows the same order over one panel each, plus `clusterNorm` in second
place: 17 panels.

Exactly **two** quartet measures are plotted, both weighted, so that the
contrast is single-axis: `wNrqs` (non-redundant quartet statements,
chance-corrected) comes **first**, ahead of `wQuartet` (naive quartet currency,
uncorrected). The unweighted `quartet` column is computed and cached but not
plotted. IQ-TREE's own `sCF` is an independent measure, not a third variant of
these.

### Fig 3 — NID vs. support scatter
Function `.NidPanel(values, name)` plots NID (x) against a support metric (y).
Points are coloured by which methods carry data for the edge (`hereFor`); the
two trend lines are GAM fits, over `common` edges and over all edges with data.
Three C-indices are chipped in each panel — no intervals, for want of space.
Layout: 6 rows × 3 cols, 17 metric panels + 1 legend slot. Output:
`Fig 3 - CID vs support.pdf` (7 × 9 in). The panels are driven by
`.nidMetrics`.

**Support chips** (`.NidChips()`, same machinery as Fig 4 — `ChipRow()`,
`ChipFill()`, `ChipInk()`, `RoundRect()`). Three chipped rows down the
top-right corner *inside* each panel;
`nidChipCex = 0.78`. Each of the three C-indices gets its own bounds on the
shared cividis ramp, all anchored on 0.5:

| Row | Subset | Domain | Observed |
|-----|--------|--------|----------|
| `C-index` | all edges with data for the measure | `cIndexRange` = 0.5–0.9 | 0.62–0.86 |
| `C_CTBI` | `common` edges only | `nidCtbiRange` = 0.5–0.95 | 0.56–0.93 |
| `C_CTBI'` | `common & !partCorrect` | `nidCtbi1Range` = 0.3–0.7 | 0.38–0.67 |

- The C-index row shares Fig 4's **domain**, but not its **subset**, so a chip
  does *not* always mean the same thing in both figures. Fig 4's `tntCF` group
  requires all five TNT statistics to be non-NA (n = 36 761); this figure's `ok`
  requires only the plotted metric. The six TNT metrics differ between the two
  by up to 0.0230 (`symFq` 0.8335 here vs 0.8565 in Fig 4) — 2.7× the widest 95%
  half-width. All ten non-TNT metrics agree to 1e-4.
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
  beneath the "Stats available" point key, which is set at `cex = 0.9` to leave
  room.
- The chips are opaque and sit over the data.

The domains above can be checked against a rendered PDF without re-running the
pipeline: cairo writes literal strings, so the printed values sit in the text
layer.

## Concordance Metrics (TreeSearch package)

| Column | Function | Notes |
|--------|----------|-------|
| `cluster` | `ClusteringConcordance(chanceCorrect = FALSE)` | Primary metric |
| `mutual` | `MutualClusteringConcordance()` | Primary metric |
| `wQuartet` | `QuartetConcordance(unit = "quartet", chanceCorrect = FALSE)` | Naive quartet currency |
| `wNrqs` | `QuartetConcordance(unit = "nrqs", chanceCorrect = TRUE)` | May be negative |
| `quartet` | as `wQuartet`, `weight = FALSE` | Cached; not plotted |
| `clusterNorm` | `ClusteringConcordance(chanceCorrect = TRUE)` | Fig 3 panel 2, "Adjusted clustering conc." |
| `phylo` | `PhylogeneticConcordance()` | Not used in study |
| `shared` | `SharedPhylogeneticConcordance()` | Not used in study |

## Key Dependencies

**R packages:** `TreeTools`, `TreeSearch`, `TreeDist`, `phangorn`, `ape`,
`survival` (C-index via `concordance`), `mgcv`, `pROC`, `viridisLite`, `cli`;
`ssh` for the HPC scripts only.

**External software:** TNT, MrBayes 3.2.7, IQ-TREE 3.0.1. R ≥ 4.1.

TreeSearch must be installed via `pak::pak("ms609/TreeSearch@cpp-search")` until
v2.0.0 reaches CRAN.


## The C-index (`SomersD` / `CIndex` in 80_byEdge.R)

Computed with **`survival::concordance()`**, whose O(n log n) balanced tree
takes 2 s at n = 354 638. Fig 3 needs 51 C-indices, six of them at full size, so
a quadratic ranker costs hours per run — do not substitute one.

- `SomersD()` takes the standard *error* of Dxy directly from
  `sqrt(fit$var)`, the infinitesimal-jackknife variance. It is already an error,
  not a standard deviation: dividing by `sqrt(n)` would understate every
  interval by a factor of √n (595× at full sample size).
- Intervals by panel group:

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

## Notes

- Concordance results are cached in `concordance/` and `entropy/` to avoid
  recomputation; delete cache files if partitions change. These two are
  distributed with the repository — a full rebuild is ~8 h. `roc/` and `nid/`
  are not, and rebuild on demand.
- `90_byChar.R` writes three per-character caches, distinguished by the
  `ConcFile()` suffix: `_chr`, `_chrQ` and `_cns`.
- `data-raw/nid/` caches each Fig 3 metric's three C-indices *and* its two GAM
  prediction curves. The guard checks `nidCacheVersion` and nothing else, so
  **bump that constant** after any change to the GAM family, `gamma`, the
  spline basis, the subsets the C-indices use, or the meaning of a metric
  column. A second run of the figure then costs seconds, not hours.
- **Never call `QuartetConcordance()` or `ClusteringConcordance()` without
  pinning `unit` and `chanceCorrect`.** TreeSearch defaults to
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
- **999, not 1 000, replicates reach the edge-wise analyses.** `gam0615` has no
  TNT output, and `80_byEdge.R` `next`s past any replicate that lacks one.
  `gam0011` has no MrBayes consensus tree, which nothing reads. A full re-run of
  the inference would fill both gaps and shift every statistic slightly.
- Site concordance factors come from the `.cf.tree` files, which ship with the
  repository. `80_byEdge.R` regenerates a missing one by calling the IQ-TREE
  executable named in the **`iqtree.exe` environment variable** — not `iqExec`
  from `_config.R`. If it is unset, `system2("")` fails quietly and the script
  then errors on the absent file.
- **Only the simulation is seeded.** `10_simulate.R` sets `set.seed(1984)`,
  IQ-TREE takes `-seed 1`, and Fig 3's scatter subsample uses `set.seed(4917)`.
  `mb-gam.nex` sets no MrBayes seed, `tnt-ew.run` sets no `rseed`, and
  `Consistency(nRelabel = 1000)` in `90_byChar.R` randomises with no seed at
  all — so re-running those steps reproduces the study's findings, not its
  digits. The distributed inference output and caches are the published record.
- **Fig 5 rests on a selected subset.** 12.41% of `_chr.txt` values (11 910 /
  96 000) are `NaN`: an invariant character gives `charInfo = charMax = 0` and
  the `return = "char"` branch has no zero guard. `na.omit()` in `90_byChar.R`
  drops them silently, and the loss is rate-correlated — 58.2% of characters in
  the slowest gamma category against 0% in the two fastest — so the leftmost box
  rests on ~42% of its characters, selected for being variable.
- The alignments contain **no missing data** and reach at most **4 states** per
  character; `nCats = 6` counts gamma rate categories, not states.
- `90_byChar.R` writes with `pdf()` where `80_byEdge.R` uses `cairo_pdf()`, so
  Fig 5 will not embed the Unicode the other figures rely on.
- `.zenodo.json` supplies the deposit's title, author, ORCID and licence. It
  carries no `related_identifiers` yet; add the article DOI as `isSupplementTo`
  once it is issued.
