# LDpred simulations

Self-contained simulation studies for design questions that the golden-file
tests cannot answer. They need only `numpy` and an installed `ldpred`.

## `lrld_policy.py` — how should long-range LD regions be handled?

LDpred models LD as a band of ±`ld_radius` variants. Inside a long-range LD
region — the MHC, a polymorphic inversion, a centromere — real LD extends far
past that band, so the model is misspecified exactly there. LDpred ships the
Price et al. (AJHG 2008) table of 34 such regions, lifted to hg38: **~168 Mb,
about 5.4% of the genome**. `ldpred gibbs` excludes them unless
`--incl-long-range-ld` is passed.

Four policies are possible:

| policy | behaviour |
|---|---|
| `include` | sample those variants normally — this is `--incl-long-range-ld` |
| `zero` | skip them in the sweep, emit weight 0 |
| `inf` | skip them in the sweep, emit their LDpred-inf weight — **what LDpred now does** |
| `drop` | remove them from the panel before LD is built |

Under either skip policy the variant is frozen at its LDpred-inf effect
*internally*, so it still absorbs LD from its neighbours. The two differ only in
what gets written out. LDpred used to emit 0, contradicting its own `--help`
text and silently dropping ~5% of the genome from the score; on the strength of
the results below it now emits the LDpred-inf weight.

The simulation sets the region's weights explicitly in every arm rather than
inheriting whatever the sampler currently emits, so `zero` stays a real
comparison rather than becoming a duplicate of `inf`.

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

| policy | loading 0.0 *(control)* | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|
| `include` | **0.435** | 0.014 | 0.029 | 0.080 |
| `zero` | 0.324 | 0.328 | 0.329 | 0.330 |
| `inf` | **0.444** | **0.374** | 0.300 | 0.288 |
| `drop` | 0.347 | 0.332 | **0.336** | **0.332** |

**Sum of squared effects** (true h² = 0.5; far above it means the sampler diverged)

| policy | 0.0 | 0.5 | 1.0 | 2.0 |
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

3. **`inf` is the only policy that is safe in both regimes.** It never diverges,
   and it pays no exclusion penalty when the region is benign (0.444 vs 0.435
   for `include`). `zero` and `drop` are safe but give up ~0.1 R² whenever the
   region turns out fine. Under strong long-range LD, `zero`/`drop` edge ahead
   of `inf` by ~0.03–0.04.

**Implication, now implemented:** LDpred emits the LDpred-inf weight for skipped
variants rather than zero, which is what its help text always claimed. A
three-way option (`include` / `inf` / `drop`) would express the trade-off better
than the current boolean, and is not implemented — under strong long-range LD
`drop` is the best of the four, by ~0.03–0.04 over `inf`.

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
