"""How should LDpred handle long-range LD regions?

LDpred models LD as a band of +/- ``ld_radius`` variants. Inside a long-range LD
region -- the MHC, a polymorphic inversion, a centromere -- real LD extends far
past that band, so the model is misspecified exactly there. LDpred ships a table
of 34 such regions (Price et al. AJHG 2008, lifted to hg38) spanning ~168 Mb,
about 5.4% of the genome, and ``ldpred gibbs`` excludes them unless
``--incl-long-range-ld`` is passed.

Four policies are possible, and the choice is empirical:

``include``
    Sample those variants normally. This is ``--incl-long-range-ld``.
``zero``
    Skip them in the sweep and emit a weight of 0. This is what the code does
    today: a skipped variant is frozen at its LDpred-inf effect internally, so it
    still absorbs LD from its neighbours, but its own posterior mean is never
    written and stays 0.
``inf``
    Skip them in the sweep but emit their LDpred-inf weight. This is what
    ``--help`` currently claims happens.
``drop``
    Remove them from the panel before LD is built, so they are absent from
    beta_hat, from every LD window, and from the score.

Run::

    python simulations/lrld_policy.py --reps 10
    python simulations/lrld_policy.py --reps 8 --loadings 0 0.5 1 2 --csv out.csv

Model. Genotypes are standardised Gaussians rather than allele counts. LDpred
standardises genotypes internally and the question here is purely about LD
structure, so this isolates the effect under test at the cost of realism about
allele frequencies. Background LD is an AR(1) process along the variant index,
which decays inside the radius. The long-range region adds a shared latent
factor, so every pair inside it is correlated regardless of distance -- the
inversion/MHC pattern. ``--loading 0`` disables the factor, which is the control
condition: a region that is on the exclusion list but behaves normally.

Choosing parameters. ``--rho`` and ``--ld-radius`` are coupled: the band has to
capture the *background* LD, or the model is misspecified everywhere and every
policy diverges, which tells you nothing about the region. Background
correlation at the band edge is ``rho ** ld_radius`` and wants to be small --
the defaults give 0.95 ** 50 = 0.08. Raising ``--rho`` or lowering
``--ld-radius`` without checking that product produces uninterpretable runs. The
printed "band captures X% of the region's LD mass" line is the diagnostic: with
``--loading 0`` it should be high (~77%), and it falls as the loading rises.

Caveats. One region, one factor, one heritability, one polygenicity. Real MHC LD
is multi-haplotype and messier. The cost of excluding scales with how much true
signal sits in the excluded regions, which ``--frac-causal-in-lrld`` controls.
Before drawing conclusions for production, run this on real genotypes with the
actual Price regions.
"""

import argparse
import csv
import sys

import numpy as np

from ldpred import LDpred_gibbs, LDpred_inf, ld

POLICIES = ('include', 'zero', 'inf', 'drop')


def simulate_panel(rng, n, m, rho, lrld_slice, lrld_loading):
    """Standardised genotypes: AR(1) background plus an optional long-range factor."""
    Z = np.empty((n, m))
    Z[:, 0] = rng.standard_normal(n)
    innovation = np.sqrt(1.0 - rho ** 2)
    for j in range(1, m):
        Z[:, j] = rho * Z[:, j - 1] + innovation * rng.standard_normal(n)
    if lrld_loading > 0:
        factor = rng.standard_normal(n)
        Z[:, lrld_slice] += lrld_loading * factor[:, None]
    Z -= Z.mean(0)
    Z /= Z.std(0)
    return Z


def build_ld(panel, ld_radius):
    """LD dict and reference matrices, in the layout the sampler expects."""
    snps = np.ascontiguousarray(panel.T, dtype=np.float32)      # (m, n)
    return ld.get_LDpred_ld_tables(snps, ld_radius=ld_radius,
                                   ld_window_size=2 * ld_radius)


def describe_ld(panel, lrld_slice, ld_radius):
    """How much of the region's LD mass does the band actually see?"""
    sub = panel[:, lrld_slice]
    R = np.corrcoef(sub, rowvar=False)
    k = R.shape[0]
    band = np.abs(np.arange(k)[:, None] - np.arange(k)[None, :]) <= ld_radius
    total = np.abs(R).sum()
    return {
        'mean_abs_r_within_band': float(np.abs(R[band]).mean()),
        'mean_abs_r_outside_band': float(np.abs(R[~band]).mean()),
        'band_share_of_ld_mass': float(np.abs(R[band]).sum() / total),
    }


def weights_for_policy(policy, beta_hat, ld_ret, ref_panel, lrld_mask,
                       m, ld_radius, h2, n_gwas, p, num_iter, burn_in):
    """LDpred weights over all m variants under one policy."""
    inf_betas = LDpred_inf.ldpred_inf(
        beta_hat, genotypes=None, reference_ld_mats=ld_ret['ref_ld_matrices'],
        h2=h2, n=n_gwas, ld_window_size=2 * ld_radius, verbose=False)

    if policy == 'drop':
        keep = ~lrld_mask
        sub_ld = build_ld(ref_panel[:, keep], ld_radius)
        sub_inf = LDpred_inf.ldpred_inf(
            beta_hat[keep], genotypes=None,
            reference_ld_mats=sub_ld['ref_ld_matrices'],
            h2=h2, n=n_gwas, ld_window_size=2 * ld_radius, verbose=False)
        res = LDpred_gibbs.ldpred_gibbs(
            beta_hat[keep], start_betas=sub_inf.copy(), h2=h2, n=n_gwas, ns=None,
            ld_radius=ld_radius, num_iter=num_iter, burn_in=burn_in, p=p,
            ld_dict=sub_ld['ld_dict'], verbose=False, print_progress=False)
        w = np.zeros(m)
        w[keep] = res['betas']
        return w

    res = LDpred_gibbs.ldpred_gibbs(
        beta_hat, start_betas=inf_betas.copy(), h2=h2, n=n_gwas, ns=None,
        ld_radius=ld_radius, num_iter=num_iter, burn_in=burn_in, p=p,
        ld_dict=ld_ret['ld_dict'],
        snp_lrld=None if policy == 'include' else lrld_mask,
        verbose=False, print_progress=False)
    w = res['betas'].copy()
    # Set the region's weights explicitly rather than relying on whatever the
    # sampler currently emits for skipped variants.  That default is the thing
    # under comparison here, so the arms must not inherit it.
    if policy == 'zero':
        w[lrld_mask] = 0.0
    elif policy == 'inf':
        w[lrld_mask] = inf_betas[lrld_mask]
    return w


def one_replicate(seed, cfg):
    rng = np.random.default_rng(seed)
    m = cfg['m']
    lrld_slice = slice(cfg['lrld_start'], cfg['lrld_start'] + cfg['lrld_width'])
    lrld_mask = np.zeros(m, dtype=bool)
    lrld_mask[lrld_slice] = True

    panel_kw = dict(rho=cfg['rho'], lrld_slice=lrld_slice,
                    lrld_loading=cfg['loading'])
    ref = simulate_panel(rng, cfg['n_ref'], m, **panel_kw)
    gwas = simulate_panel(rng, cfg['n_gwas'], m, **panel_kw)
    test = simulate_panel(rng, cfg['n_test'], m, **panel_kw)

    # Causal variants, deliberately enriched inside the region so the policies
    # are not separated trivially by "there was no signal there anyway".
    n_causal = cfg['n_causal']
    n_in = int(round(n_causal * cfg['frac_causal_in_lrld']))
    idx_in = rng.choice(np.flatnonzero(lrld_mask), n_in, replace=False)
    idx_out = rng.choice(np.flatnonzero(~lrld_mask), n_causal - n_in, replace=False)
    beta_true = np.zeros(m)
    beta_true[np.concatenate([idx_in, idx_out])] = (
        rng.standard_normal(n_causal) * np.sqrt(cfg['h2'] / n_causal))

    def phenotype(X):
        g = X @ beta_true
        g *= np.sqrt(cfg['h2']) / g.std()
        return g + rng.standard_normal(len(g)) * np.sqrt(1.0 - cfg['h2'])

    y_gwas = phenotype(gwas)
    y_test = phenotype(test)
    y_test = (y_test - y_test.mean()) / y_test.std()

    yg = (y_gwas - y_gwas.mean()) / y_gwas.std()
    beta_hat = (gwas.T @ yg) / cfg['n_gwas']        # marginal, standardised scale

    ld_ret = build_ld(ref, cfg['ld_radius'])

    row = {'seed': seed}
    row.update(describe_ld(ref, lrld_slice, cfg['ld_radius']))
    for policy in POLICIES:
        w = weights_for_policy(
            policy, beta_hat, ld_ret, ref, lrld_mask, m, cfg['ld_radius'],
            cfg['h2'], cfg['n_gwas'], n_causal / m, cfg['num_iter'], cfg['burn_in'])
        prs = test @ w
        row[policy + '_r2'] = (np.corrcoef(prs, y_test)[0, 1] ** 2
                               if prs.std() > 0 else 0.0)
        # Against the true h2, this is LDpred's own convergence diagnostic:
        # a diverged sampler shows a sum of squared effects far above h2.
        row[policy + '_sum_sq'] = float(np.sum(w ** 2))
    return row


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.split('\n\n')[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--reps', type=int, default=10, help='replicates per loading')
    p.add_argument('--loadings', type=float, nargs='+', default=[0.0, 0.5, 1.0, 2.0],
                   help='long-range factor loadings; 0 is the control')
    p.add_argument('--m', type=int, default=4000, help='number of variants')
    p.add_argument('--ld-radius', type=int, default=50)
    p.add_argument('--lrld-start', type=int, default=1500)
    p.add_argument('--lrld-width', type=int, default=600,
                   help='width of the long-range region, in variants')
    p.add_argument('--rho', type=float, default=0.95, help='AR(1) background LD')
    p.add_argument('--n-ref', type=int, default=2000)
    p.add_argument('--n-gwas', type=int, default=20000)
    p.add_argument('--n-test', type=int, default=3000)
    p.add_argument('--h2', type=float, default=0.5)
    p.add_argument('--n-causal', type=int, default=40)
    p.add_argument('--frac-causal-in-lrld', type=float, default=0.25)
    p.add_argument('--num-iter', type=int, default=60)
    p.add_argument('--burn-in', type=int, default=10)
    p.add_argument('--csv', default=None, help='write per-replicate rows here')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    all_rows = []
    for loading in args.loadings:
        cfg = dict(vars(args))
        cfg['loading'] = loading
        cfg['ld_radius'] = args.ld_radius
        rows = [one_replicate(seed, cfg) for seed in range(args.reps)]
        for r in rows:
            r['loading'] = loading
        all_rows.extend(rows)

        band_share = np.mean([r['band_share_of_ld_mass'] for r in rows])
        print('\n  long-range factor loading = %g   (%d replicates)' % (loading, args.reps))
        print('  band captures %.1f%% of the region\'s LD mass' % (100 * band_share))
        print('  %-9s %14s %9s %14s' % ('policy', 'mean test R2', 'se', 'mean sum b^2'))
        base = np.array([r['include_r2'] for r in rows])
        for policy in POLICIES:
            v = np.array([r[policy + '_r2'] for r in rows])
            ssq = np.mean([r[policy + '_sum_sq'] for r in rows])
            se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else float('nan')
            delta = ''
            if policy != 'include' and len(v) > 1:
                d = v - base
                delta = '   delta vs include = %+.4f (se %.4f)' % (
                    d.mean(), d.std(ddof=1) / np.sqrt(len(d)))
            print('  %-9s %14.4f %9.4f %14.2f%s' % (policy, v.mean(), se, ssq, delta))
        print('  (true h2 = %g; sum b^2 far above it means the sampler diverged)' % args.h2)

    if args.csv:
        fields = sorted({k for r in all_rows for k in r})
        with open(args.csv, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(all_rows)
        print('\n  wrote %d rows to %s' % (len(all_rows), args.csv))
    return 0


if __name__ == '__main__':
    sys.exit(main())
