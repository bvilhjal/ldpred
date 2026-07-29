# LDpred simulations

Self-contained simulation studies for design questions that the golden-file
tests cannot answer. They need only `numpy` and an installed `ldpred`.

## `lrld_policy.py` — how should long-range LD regions be handled?

LDpred models LD as a band of ±`ld_radius` variants. Inside a long-range LD
region — the MHC, a polymorphic inversion, a centromere — real LD extends far
past that band, so the model is misspecified exactly there. LDpred ships the
Price et al. (AJHG 2008) table of 34 such regions, lifted to hg38: **~168 Mb,
about 5.4% of the genome**. `ldpred gibbs` holds them out of the sampler unless
told otherwise, via `--long-range-ld {include,inf,exclude}`.

Five arms are compared here; three of them are reachable from the command line.

| arm | behaviour | exposed as |
|---|---|---|
| `include` | sample those variants normally | `--long-range-ld include` |
| `zero` | skip them, keep them in the model at their LDpred-inf effect, emit weight 0 | — (the old behaviour) |
| `inf` | skip them, keep them in the model, emit their LDpred-inf weight | `--long-range-ld inf` **(default)** |
| `exclude` | skip them with a zero effect, so they leave the model, emit weight 0 | `--long-range-ld exclude` |
| `drop` | remove them from the panel before LD is built | — (coordination-step only) |

A held-out variant stays frozen at whatever value it starts with, so that
starting value decides two things at once: what it contributes to its
neighbours' residuals, and what weight it is reported with. `zero` and `inf`
keep it at its LDpred-inf effect, so it still absorbs LD from its neighbours,
and differ only in what gets written out. `exclude` starts it at zero, so it
contributes nothing and leaves the model altogether.

`zero` was LDpred's behaviour before this work: it contradicted the `--help`
text and silently dropped ~5% of the genome from the score. It is kept here as
a comparison arm but is no longer reachable.

The simulation sets each arm's region weights explicitly rather than inheriting
whatever the sampler currently emits, so the arms stay independent of the
default.

### Running

```bash
python simulations/lrld_policy.py --reps 10
python simulations/lrld_policy.py --reps 10 --loadings 0 0.5 1 2 --csv results.csv
```

`--loading` controls how much long-range LD the region has. **`--loading 0` is
the control**: a region that sits on the exclusion list but behaves normally.

### Results

10 replicates per cell, defaults otherwise. Full per-replicate output in
`lrld_policy_results.csv`.

**Mean test R²**

| arm | loading 0.0 *(control)* | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|
| `include` | **0.435** | 0.014 | 0.029 | 0.080 |
| `zero` | 0.324 | 0.328 | 0.329 | 0.330 |
| `inf` | **0.444** | **0.374** | 0.300 | 0.288 |
| `exclude` | 0.326 | 0.333 | **0.334** | **0.333** |
| `drop` | 0.347 | 0.332 | **0.336** | **0.332** |

**Sum of squared effects** (true h² = 0.5; far above it means the sampler diverged)

| arm | 0.0 | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|
| `include` | 2.2 | **429.9** | **244.3** | **126.9** |
| others | ~0.4–1.5 | ~0.5 | ~0.5 | ~0.5–0.8 |

### What this shows

1. **`include` does not merely do worse — the sampler diverges.** Sum of squared
   effects reaches 250–860× the true heritability, concentrated in the region.
   The band captures only ~21% of the region's LD mass at loading 1.0, so the
   local system is near-singular. LDpred's existing warning — *"Sum of squared
   updated effects estimates seems too large… suggests that the Gibbs sampler
   did not converge"* — is the signal to watch for on real data.

2. **Excluding is not free.** In the control column, `zero` costs **−0.111 R²**
   and `drop` **−0.087** against including. That matters because the Price list
   is broad and derived from 2008 European data; plenty of it may be benign for
   a given trait or ancestry.

3. **`inf` is the only arm that is safe in both regimes.** It never diverges,
   and it pays no exclusion penalty when the region is benign (0.444 vs 0.435
   for `include`). `zero`, `exclude` and `drop` are all safe but give up
   ~0.09–0.12 R² whenever the region turns out fine. Under strong long-range
   LD they edge ahead of `inf` by ~0.03–0.04.

4. **`exclude` is a faithful stand-in for `drop`.** Wherever there is real
   long-range LD it matches the full panel rebuild to within noise (0.334 vs
   0.336 at loading 1.0; 0.333 vs 0.332 at 2.0). It is cheaper because it needs
   no re-derivation of LD, which is why it is the one exposed on the command
   line. It falls short of `drop` only in the control column (0.326 vs 0.347),
   which is the regime where `inf` is the right choice anyway.

**Implemented as a three-way option**, `--long-range-ld {include,inf,exclude}`,
defaulting to `inf`:

- `inf` when you do not know whether a listed region really carries long-range
  LD in your data — it is the only arm that is never bad.
- `exclude` when the regions genuinely do carry strong long-range LD — it is
  the best of the reachable arms at loading >= 0.5.
- `include` only to reproduce the old (pre-fix) behaviour; expect divergence.

`drop` itself is not exposed: rebuilding every retained variant's LD window on
the reduced panel is a coordination-step operation, not something the sampler
can do, and `exclude` captures nearly all of its benefit.

### Model and caveats

Genotypes are standardised Gaussians, not allele counts — LDpred standardises
internally and the question is purely about LD structure, so this isolates the
effect at the cost of realism about allele frequencies. Background LD is an
AR(1) process along the variant index; the long-range region adds a shared
latent factor, so every pair inside it is correlated regardless of distance.

- One region, one factor, one heritability, one polygenicity. Real MHC LD is
  multi-haplotype and messier.
- The cost of excluding scales with how much true signal sits in the excluded
  regions (`--frac-causal-in-lrld`, default 0.25 — deliberate enrichment).
- `--rho` and `--ld-radius` are coupled. The band must capture the *background*
  LD or everything diverges and the run says nothing about the region.
  Background correlation at the band edge is `rho ** ld_radius`; the defaults
  give 0.95 ** 50 = 0.08. The printed "band captures X% of the region's LD mass"
  line is the diagnostic.

**Before changing production behaviour, this should be repeated on real
genotypes with the actual Price regions.**
