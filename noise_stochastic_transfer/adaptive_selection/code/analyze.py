"""Pre-specified analysis of results/raw.json (PROTOCOL.md section 3)."""
import os, json, math, csv
import numpy as np
from scipy.stats import beta, norm
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
D = json.load(open(os.path.join(RES, 'raw.json')))
ALPHA = D['alpha']
CELLS = D['configs']
RULES = ['S', 'I', 'C', 'L', 'C2']
NAME = {'S': 'S (selection)', 'I': 'I (agreement)', 'C': 'C (cost)', 'L': 'L (linear)', 'C2': 'C½ (cost, α/2)'}
XPOS = {'m2': -2.0, 'm0.5': -0.5, 'm0.1': -0.1, 'between': None, 'p0.1': 0.1, 'p0.5': 0.5, 'p2': 2.0}
POSLAB = {'m2': 'λc−2w', 'm0.5': 'λc−0.5w', 'm0.1': 'λc−0.1w', 'between': '(λ*+λc)/2', 'p0.1': 'λc+0.1w', 'p0.5': 'λc+0.5w', 'p2': 'λc+2w'}
ORDER_POS = ['m2', 'm0.5', 'm0.1', 'between', 'p0.1', 'p0.5', 'p2']
NCELL = len(CELLS)
ZB = norm.ppf(1 - 0.05 / NCELL / 2)


def wilson(k, n, zq=1.959964):
    if n == 0:
        return (float('nan'),) * 3
    p = k / n; den = 1 + zq ** 2 / n
    c = (p + zq ** 2 / (2 * n)) / den; h = zq * math.sqrt(p * (1 - p) / n + zq ** 2 / (4 * n ** 2)) / den
    return p, max(0.0, c - h), min(1.0, c + h)


def cp(k, n):
    if n == 0:
        return (float('nan'),) * 3
    lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
    return k / n, lo, hi


def dec(r, k):
    return r['d' + k]


rows = []
violations_I_sub_S = 0
for c in sorted(CELLS, key=lambda c: (['pair', 'sym3', 'asym3'].index(c['cluster']), -c['z'], c['nu'], ORDER_POS.index(c['pos']))):
    R, tr = c['records'], c['truth']; n = len(R)
    row = dict(cluster=c['cluster'], z=c['z'], nu=c['nu'], pos=c['pos'], truth=tr, n=n,
               x_over_w=(c['lam'] - c['lam_c']) / c['width'], lamstar_minus_lamc_over_w=(c['lam_star'] - c['lam_c']) / c['width'])
    for k in RULES:
        corr = sum(dec(r, k) == tr for r in R); wr = sum(dec(r, k) not in ('inc', tr) for r in R)
        row[k + '_correct'], row[k + '_wrong'], row[k + '_abstain'] = corr, wr, n - corr - wr
        row[k + '_wrong_rate'] = wilson(wr, n)
        row[k + '_flag'] = wilson(wr, n)[1] > ALPHA
        row[k + '_flag_bonf'] = wilson(wr, n, ZB)[1] > ALPHA
        row[k + '_wrong_given_decl'] = cp(wr, corr + wr)
    row['sel_C'] = sum(r['sel'] == 'C' for r in R); row['sel_L'] = n - row['sel_C']
    row['sel_exceptional'] = sum(not (r['A_hat'] > 0 and r['m0_hat'] > 0) for r in R)
    # paired S vs I
    for r in R:
        if r['dI'] != 'inc' and r['dS'] != r['dI']:
            violations_I_sub_S += 1
    row['S_restored_correct_vs_I'] = sum(r['dS'] == tr and r['dI'] == 'inc' for r in R)
    row['S_added_wrong_vs_I'] = sum(r['dS'] not in ('inc', tr) and r['dI'] == 'inc' for r in R)
    # paired S vs C
    row['correct_S_only_vs_C'] = sum(r['dS'] == tr and r['dC'] != tr for r in R)
    row['correct_C_only_vs_S'] = sum(r['dC'] == tr and r['dS'] != tr for r in R)
    row['wrong_S_only_vs_C'] = sum(r['dS'] not in ('inc', tr) and not (r['dC'] not in ('inc', tr)) for r in R)
    row['wrong_C_only_vs_S'] = sum(r['dC'] not in ('inc', tr) and not (r['dS'] not in ('inc', tr)) for r in R)
    # fitted optimizer, separate
    row['fit_out'] = sum(r['fit'] == 'out' for r in R)
    row['fit_wrong'] = sum(r['fit'] in ('A', 'B') and r['fit'] != tr for r in R)
    for k in RULES:
        agree = [r for r in R if dec(r, k) != 'inc' and dec(r, k) == r['fit']]
        row[k + '_eq_fit'] = len(agree); row[k + '_eq_fit_fitwrong'] = sum(r['fit'] != tr for r in agree)
    # diagnostics
    ok = [r for r in R if r['A_hat'] > 0 and r['m0_hat'] > 0]
    row['Ahat_rel_diff_max'] = max(abs(r['A_hat'] - r['A_hat_L']) / max(abs(r['A_hat_L']), 1e-300) for r in R)
    row['nuhat_over_nu_median'] = float(np.median([r['nu_hat'] / c['nu'] for r in ok])) if ok else float('nan')
    row['switch_nominal'] = (c['A'] / c['mu0']) ** D['gamma'][str(1 if c['cluster'] == 'pair' else 2)]
    rows.append(row)

fmt = lambda t: ('%.3f [%.3f, %.3f]' % t) if t[0] == t[0] else 'undefined'
with open(os.path.join(RES, 'per_cell.csv'), 'w', newline='') as f:
    w = csv.writer(f); keys = list(rows[0].keys()); w.writerow(keys)
    for r in rows:
        w.writerow([fmt(r[k]) if isinstance(r[k], tuple) else r[k] for k in keys])

CN = {'pair': 'pair', 'sym3': 'sym3', 'asym3': 'asym3'}
with open(os.path.join(RES, 'per_cell.md'), 'w') as f:
    f.write('# Per-cell results (fresh draws, PROTOCOL.md)\n\n')
    f.write('Each entry: correct/wrong/abstain out of 1000. "wrong|decl" = wrong declarations / declarations (Clopper–Pearson 95%). '
            '"uncond" = wrong/1000 (Wilson 95%); † = Wilson lower bound above α; ‡ = also above α at the Bonferroni level. '
            '"sel C" = draws where S selected the cost test.\n\n')
    f.write('| case | ν | position | truth | sel C | S ✓/✗/– | S uncond | S wrong\\|decl | I ✓/✗/– | I uncond | C ✓/✗/– | C uncond | L ✓/✗/– | L uncond | C½ ✓/✗/– |\n')
    f.write('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        def cell(k):
            return '%d/%d/%d' % (r[k + '_correct'], r[k + '_wrong'], r[k + '_abstain'])
        def un(k):
            return '%.3f%s' % (r[k + '_wrong_rate'][0], '‡' if r[k + '_flag_bonf'] else ('†' if r[k + '_flag'] else ''))
        wd = r['S_wrong_given_decl']
        wds = ('%d/%d = %.2f [%.2f, %.2f]' % (r['S_wrong'], r['S_wrong'] + r['S_correct'], *wd)) if wd[0] == wd[0] else '0/0 undefined'
        f.write('| %s z=%g | %g | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n' % (
            CN[r['cluster']], r['z'], r['nu'], POSLAB[r['pos']], r['truth'], r['sel_C'], cell('S'), un('S'), wds,
            cell('I'), un('I'), cell('C'), un('C'), cell('L'), un('L'), cell('C2')))

with open(os.path.join(RES, 'paired.md'), 'w') as f:
    f.write('# Paired comparisons on the same draws\n\n')
    f.write('I ⊂ S draw by draw (Proposition 3(a)): violations = %d.\n\n' % violations_I_sub_S)
    f.write('| case | ν | position | truth | S restores correct (vs I) | S adds wrong (vs I) | correct S only / C only | wrong S only / C only |\n|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        f.write('| %s z=%g | %g | %s | %s | %d %s | %d %s | %d / %d | %d / %d |\n' % (
            CN[r['cluster']], r['z'], r['nu'], POSLAB[r['pos']], r['truth'],
            r['S_restored_correct_vs_I'], '[%.3f, %.3f]' % wilson(r['S_restored_correct_vs_I'], r['n'])[1:],
            r['S_added_wrong_vs_I'], '[%.3f, %.3f]' % wilson(r['S_added_wrong_vs_I'], r['n'])[1:],
            r['correct_S_only_vs_C'], r['correct_C_only_vs_S'], r['wrong_S_only_vs_C'], r['wrong_C_only_vs_S']))

with open(os.path.join(RES, 'fit_agreement.md'), 'w') as f:
    f.write('# Agreement with the fitted optimizer (separate from branch declarations)\n\n')
    f.write('Fits are numerical multistart optima, not certified global optima. "k=fit (fit wrong)" counts draws where rule k declares '
            'and the complex global fit lies in the declared arc, and among them those where the fitted branch is not L(x).\n\n')
    f.write('| case | ν | position | truth | fit outside arcs | fit wrong | S=fit (fit wrong) | I=fit (fit wrong) | C=fit (fit wrong) |\n|---|---|---|---|---|---|---|---|---|\n')
    for r in rows:
        f.write('| %s z=%g | %g | %s | %s | %d | %d | %d (%d) | %d (%d) | %d (%d) |\n' % (
            CN[r['cluster']], r['z'], r['nu'], POSLAB[r['pos']], r['truth'], r['fit_out'], r['fit_wrong'],
            r['S_eq_fit'], r['S_eq_fit_fitwrong'], r['I_eq_fit'], r['I_eq_fit_fitwrong'], r['C_eq_fit'], r['C_eq_fit_fitwrong']))

summ = dict(cells=len(rows), draws=sum(r['n'] for r in rows), I_sub_S_violations=violations_I_sub_S)
for k in RULES:
    summ[k] = dict(flagged=[(r['cluster'], r['z'], r['nu'], r['pos'], round(r[k + '_wrong_rate'][0], 3)) for r in rows if r[k + '_flag']],
                   flagged_bonf=[(r['cluster'], r['z'], r['nu'], r['pos'], round(r[k + '_wrong_rate'][0], 3)) for r in rows if r[k + '_flag_bonf']],
                   max_wrong=max(r[k + '_wrong_rate'][0] for r in rows),
                   pooled_correct=sum(r[k + '_correct'] for r in rows), pooled_wrong=sum(r[k + '_wrong'] for r in rows))
summ['restored_correct_vs_I_total'] = sum(r['S_restored_correct_vs_I'] for r in rows)
summ['added_wrong_vs_I_total'] = sum(r['S_added_wrong_vs_I'] for r in rows)
summ['Ahat_selector_vs_L_max_rel_diff'] = max(r['Ahat_rel_diff_max'] for r in rows)
summ['sel_exceptional_total'] = sum(r['sel_exceptional'] for r in rows)
json.dump(summ, open(os.path.join(RES, 'summary.json'), 'w'), indent=1)
json.dump(rows, open(os.path.join(RES, 'rows.json'), 'w'), default=lambda o: o.item() if hasattr(o, 'item') else str(o))

# Figure: correct and wrong declaration rates by position, for each case and noise level.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 7, 'font.family': 'serif', 'axes.linewidth': 0.5, 'legend.frameon': False})
COL = {'S': '#2a78d6', 'I': '#eb6834', 'C': '#1baf7a', 'L': '#9b59b6', 'C2': '#7f7f7f'}
MK = {'S': 'o', 'I': 's', 'C': '^', 'L': 'D', 'C2': 'v'}
cases = [('pair', 0.3), ('sym3', 0.3), ('asym3', 0.3), ('asym3', 0.15)]
nus = sorted({r['nu'] for r in rows})
for metric, fn, ylab in (('correct', 'fig_correct.pdf', 'correct-declaration rate'), ('wrong', 'fig_wrong.pdf', 'wrong-declaration rate')):
    fig, axs = plt.subplots(len(cases), len(nus), figsize=(11, 6.2), sharex=True, sharey=True)
    for i, (cl, z) in enumerate(cases):
        for j, nu in enumerate(nus):
            ax = axs[i, j]
            sub = {r['pos']: r for r in rows if r['cluster'] == cl and r['z'] == z and r['nu'] == nu}
            xs = list(range(len(ORDER_POS)))
            for k in RULES:
                ys = [sub[p][k + '_' + metric] / sub[p]['n'] for p in ORDER_POS]
                ax.plot(xs, ys, color=COL[k], marker=MK[k], ms=2.5, lw=0.9, label=NAME[k])
            if metric == 'wrong':
                ax.axhline(ALPHA, color='k', lw=0.5, ls=':')
            ax.axvline(3, color='#cccccc', lw=0.5)
            if i == 0:
                ax.set_title('ν = %g' % nu)
            if j == 0:
                ax.set_ylabel('%s z=%g\n%s' % (cl, z, ylab if i == 1 else ''))
            ax.set_xticks(xs); ax.set_xticklabels(['−2', '−.5', '−.1', 'btw', '+.1', '+.5', '+2'], fontsize=5.5)
    axs[0, 0].legend(fontsize=5.5, loc='upper left')
    fig.supxlabel('position (λ − λc)/w; "btw" = (λ*+λc)/2', fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, fn)); fig.savefig(os.path.join(RES, fn.replace('.pdf', '.png')), dpi=160)
print(json.dumps({k: (v if k not in RULES else {kk: vv for kk, vv in v.items() if kk != 'flagged'}) for k, v in summ.items()}, indent=1)[:3000])
