"""Figures and tables of the regime-selection manuscript, from the stored records of the
pre-specified study (../../../adaptive_selection/results/raw.json). No simulation is run."""
import os, json, math, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import beta, norm

HERE = os.path.dirname(os.path.abspath(__file__))
AS = os.path.join(HERE, '..', 'noise_stochastic_transfer', 'adaptive_selection')
D = json.load(open(os.path.join(AS, 'results', 'raw.json')))
ALPHA = D['alpha']; CELLS = D['configs']
RULES = ['S', 'I', 'C', 'L', 'C2']
LAB = {'S': 'selection rule S', 'I': 'agreement rule I', 'C': 'cost test C', 'L': 'linear test L', 'C2': r'cost test at $\alpha/2$'}
POS = ['m2', 'm0.5', 'm0.1', 'between', 'p0.1', 'p0.5', 'p2']
PTEX = {'m2': r'$\lambda_c-2w$', 'm0.5': r'$\lambda_c-0.5w$', 'm0.1': r'$\lambda_c-0.1w$', 'between': r'$(\lambda_*+\lambda_c)/2$',
        'p0.1': r'$\lambda_c+0.1w$', 'p0.5': r'$\lambda_c+0.5w$', 'p2': r'$\lambda_c+2w$'}
CASE = [('pair', 0.3), ('sym3', 0.3), ('asym3', 0.3), ('asym3', 0.15)]
CTEX = {('pair', 0.3): 'pair, $z=0.3$', ('sym3', 0.3): 'symmetric triple, $z=0.3$', ('asym3', 0.3): 'unequal triple, $z=0.3$', ('asym3', 0.15): 'unequal triple, $z=0.15$'}
NCELL = len(CELLS); ZB = norm.ppf(1 - 0.05 / NCELL / 2)


def wilson(k, n, zq=1.959964):
    p = k / n; den = 1 + zq ** 2 / n
    c = (p + zq ** 2 / (2 * n)) / den; h = zq * math.sqrt(p * (1 - p) / n + zq ** 2 / (4 * n ** 2)) / den
    return p, max(0.0, c - h), min(1.0, c + h)


def cp(k, n):
    if n == 0:
        return None
    return k / n, (0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1)), (1.0 if k == n else beta.ppf(0.975, k + 1, n - k))


def cell(cl, z, nu, pos):
    return next(c for c in CELLS if c['cluster'] == cl and c['z'] == z and c['nu'] == nu and c['pos'] == pos)


def counts(c, k):
    tr = c['truth']; R = c['records']
    cor = sum(r['d' + k] == tr for r in R); wr = sum(r['d' + k] not in ('inc', tr) for r in R)
    return cor, wr, len(R) - cor - wr


def selC(c):
    return sum(r['sel'] == 'C' for r in c['records'])


f3 = lambda x: ('%.3f' % x)
# ---------------- Figure 1: symmetric triple, three noise regimes ----------------
plt.rcParams.update({'font.size': 8.5, 'font.family': 'serif', 'axes.linewidth': 0.6, 'legend.frameon': False,
                     'axes.edgecolor': '#52514e', 'xtick.color': '#52514e', 'ytick.color': '#52514e'})
COL = {'S': '#2a78d6', 'I': '#eb6834', 'C': '#1baf7a', 'L': '#9b59b6', 'C2': '#8a8a8a'}
MK = {'S': 'o', 'I': 's', 'C': '^', 'L': 'D', 'C2': 'v'}
LS = {'S': '-', 'I': '--', 'C': '--', 'L': ':', 'C2': ':'}
NUF = [0.005, 0.5, 1.0]
FPOS = ['m2', 'm0.5', 'm0.1', 'p0.1', 'p0.5', 'p2']
fig, axs = plt.subplots(2, 3, figsize=(6.6, 4.0), sharex=True)
xs = np.arange(len(FPOS))
for j, nu in enumerate(NUF):
    for k in RULES:
        cs = [cell('sym3', 0.3, nu, p) for p in FPOS]
        cr = [counts(c, k)[0] / 1000 for c in cs]; wr = [counts(c, k)[1] / 1000 for c in cs]
        axs[0, j].plot(xs, cr, color=COL[k], marker=MK[k], ms=3.5, lw=1.3 if k == 'S' else 1.0, ls=LS[k], label=LAB[k])
        if k == 'S':
            wi = [wilson(counts(c, k)[1], 1000) for c in cs]
            axs[1, j].errorbar(xs, wr, yerr=[[max(0, w[0] - w[1]) for w in wi], [max(0, w[2] - w[0]) for w in wi]], color=COL[k], marker=MK[k],
                               ms=3.5, lw=1.3, ls=LS[k], elinewidth=0.7, capsize=1.5)
        else:
            axs[1, j].plot(xs, wr, color=COL[k], marker=MK[k], ms=3, lw=1.0, ls=LS[k])
    axs[1, j].axhline(ALPHA, color='k', lw=0.6, ls='-')
    axs[0, j].set_title(r'$\nu=%g$ (C selected in %d%%)' % (nu, round(100 * np.mean([selC(cell('sym3', 0.3, nu, p)) for p in FPOS]) / 1000)), fontsize=8.5)
    axs[1, j].set_ylim(0, 0.12)
    axs[0, j].set_ylim(0, 1.02)
    for a in axs[:, j]:
        a.axvline(2.5, color='#d0d0d0', lw=0.6); a.grid(True, color='#ececec', lw=0.5); a.set_axisbelow(True)
    axs[1, j].set_xticks(xs); axs[1, j].set_xticklabels(['$-2$', '$-.5$', '$-.1$', '$+.1$', '$+.5$', '$+2$'], fontsize=7.5)
axs[0, 0].set_ylabel('correct declarations'); axs[1, 0].set_ylabel('wrong declarations')
for j in (1, 2):
    axs[0, j].set_yticklabels([]); axs[1, j].set_yticklabels([])
h, l = axs[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=5, fontsize=7.2, handlelength=2.2, columnspacing=1.0, bbox_to_anchor=(0.5, -0.005))
axs[1, 1].set_xlabel(r'amplitude position $(\lambda-\lambda_c)/w$', fontsize=8)
fig.tight_layout(rect=(0, 0.06, 1, 1), h_pad=0.5, w_pad=0.6)
fig.savefig(os.path.join(HERE, 'fig_sym3.pdf')); fig.savefig(os.path.join(HERE, 'fig_sym3.png'), dpi=200)

# ---------------- Figure 2: selection frequency against noise level ----------------
NUS = sorted({c['nu'] for c in CELLS})
fig, ax = plt.subplots(figsize=(4.0, 2.4))
CC = {('pair', 0.3): '#2a78d6', ('sym3', 0.3): '#eb6834', ('asym3', 0.3): '#1baf7a', ('asym3', 0.15): '#9b59b6'}
for (cl, z) in CASE:
    fr = [np.mean([selC(cell(cl, z, nu, p)) for p in POS]) / 1000 for nu in NUS]
    ax.plot(range(len(NUS)), fr, marker='o', ms=3.5, lw=1.2, color=CC[(cl, z)], label=CTEX[(cl, z)].replace('$', ''))
ax.set_xticks(range(len(NUS))); ax.set_xticklabels([('%g' % n) for n in NUS], fontsize=7.5)
ax.set_xlabel(r'noise level $\nu$ (design levels)'); ax.set_ylabel('fraction selecting C'); ax.set_ylim(-0.02, 1.02)
ax.grid(True, color='#ececec', lw=0.5); ax.legend(fontsize=7, loc='lower left')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_selection.pdf')); fig.savefig(os.path.join(HERE, 'fig_selection.png'), dpi=200)

# ---------------- Table 1: calibration summary ----------------
rows = []
for k in RULES:
    fl = sum(wilson(counts(c, k)[1], 1000)[1] > ALPHA for c in CELLS)
    fb = sum(wilson(counts(c, k)[1], 1000, ZB)[1] > ALPHA for c in CELLS)
    worst = max(CELLS, key=lambda c: counts(c, k)[1]); wc = counts(worst, k)[1]; wi = wilson(wc, 1000)
    pc = sum(counts(c, k)[0] for c in CELLS); pw = sum(counts(c, k)[1] for c in CELLS)
    rows.append((LAB[k], fl, fb, '%s [%.4f, %.4f]' % (f3(wi[0]), wi[1], wi[2]), pc, pw))
with open(os.path.join(HERE, 'table1.tex'), 'w') as g:
    g.write('\\begin{tabular}{@{}lccccc@{}}\n\\toprule\n & & & largest & \\multicolumn{2}{c}{pooled}\\\\\n\\cmidrule(lr){5-6}\n rule & flagged & Bonf. & wrong rate [95\\%] & correct & wrong\\\\\n\\midrule\n')
    for r in rows:
        g.write('%s & %d & %d & %s & %s & %s\\\\\n' % (r[0], r[1], r[2], r[3], format(r[4], ','), format(r[5], ',')))
    g.write('\\bottomrule\n\\end{tabular}\n')

# ---------------- Table 2: focus cells ----------------
FOC = [('sym3', 0.3, 0.005, p) for p in ('m2', 'm0.5', 'm0.1', 'between')] + [('sym3', 0.3, 0.02, p) for p in ('m2', 'm0.1')] + \
      [('asym3', 0.15, 0.005, 'p0.1')]
with open(os.path.join(HERE, 'table2.tex'), 'w') as g:
    g.write('\\begin{tabular}{@{}lccccc@{}}\n\\toprule\n position (truth) & I \\checkmark/$\\times$/-- & S \\checkmark/$\\times$/-- & I wrong$\\mid$decl. & S wrong$\\mid$decl. & restored\\\\\n\\midrule\n')
    last = None
    for (cl, z, nu, p) in FOC:
        c = cell(cl, z, nu, p); a = counts(c, 'I'); b = counts(c, 'S')
        ci = '%d/%d' % (a[1], a[0] + a[1]) if a[0] + a[1] else '0/0'
        cs = '%d/%d' % (b[1], b[0] + b[1]) if b[0] + b[1] else '0/0'
        rest = sum(r['dS'] == c['truth'] and r['dI'] == 'inc' for r in c['records'])
        if (cl, z, nu) != last:
            g.write('\\multicolumn{6}{@{}l}{\\emph{%s, $\\nu=%g$}}\\\\\n' % (CTEX[(cl, z)], nu))
        last = (cl, z, nu)
        g.write('%s ($%s$) & %d/%d/%d & %d/%d/%d & %s & %s & %d\\\\\n' % (PTEX[p], c['truth'], *a, *b, ci, cs, rest))
    g.write('\\bottomrule\n\\end{tabular}\n')

# ---------------- Table 3: transition cells (nu = 0.5), near the tie ----------------
with open(os.path.join(HERE, 'table3.tex'), 'w') as g:
    g.write('\\begin{tabular}{@{}llccccc@{}}\n\\toprule\n case & position (truth) & C selected & S \\checkmark/$\\times$ & C \\checkmark/$\\times$ & L \\checkmark/$\\times$ & I \\checkmark/$\\times$\\\\\n\\midrule\n')
    for (cl, z) in CASE:
        for i, p in enumerate(('m0.1', 'between', 'p0.1')):
            c = cell(cl, z, 0.5, p)
            g.write('%s & %s ($%s$) & %d & %d/%d & %d/%d & %d/%d & %d/%d\\\\\n' % (CTEX[(cl, z)] if i == 0 else '', PTEX[p], c['truth'], selC(c),
                    *counts(c, 'S')[:2], *counts(c, 'C')[:2], *counts(c, 'L')[:2], *counts(c, 'I')[:2]))
    g.write('\\bottomrule\n\\end{tabular}\n')

# ---------------- Table 4: pooled by noise level (secondary) ----------------
with open(os.path.join(HERE, 'table4.tex'), 'w') as g:
    g.write('\\begin{tabular}{@{}lc' + 'c' * 5 + '@{}}\n\\toprule\n $\\nu$ & C selected & S & I & C & L & C at $\\alpha/2$\\\\\n\\midrule\n')
    for nu in NUS:
        cs = [c for c in CELLS if c['nu'] == nu]
        sc = sum(selC(c) for c in cs) / (1000 * len(cs))
        ent = []
        for k in RULES:
            ent.append('%s/%s' % (format(sum(counts(c, k)[0] for c in cs), ','), format(sum(counts(c, k)[1] for c in cs), ',')))
        g.write('%g & %.0f\\%% & %s\\\\\n' % (nu, 100 * sc, ' & '.join(ent)))
    g.write('\\bottomrule\n\\end{tabular}\n')

# ---------------- Supplement: per-cell table (all 196 cells) ----------------
CN = {'pair': 'pair', 'sym3': 'sym.', 'asym3': 'uneq.'}
with open(os.path.join(HERE, 'tableS_cells.tex'), 'w') as g:
    head = ('case & $\\nu$ & position & truth & sel.\\ C & S & I & C & L & C$_{\\alpha/2}$ & S wrong$\\mid$decl.\\\\\n')
    g.write('\\begin{longtable}{@{}lrllrlllllr@{}}\n\\caption{Per-cell counts (correct/wrong/abstain out of 1000 draws) for every rule, the number of draws in which S selected the cost test, and the wrong declarations of S among its declarations. A dagger marks a Wilson 95\\% lower bound above $\\alpha$.}\\label{tab:cells}\\\\\n\\toprule\n' + head + '\\midrule\n\\endfirsthead\n\\toprule\n' + head + '\\midrule\n\\endhead\n')
    for (cl, z) in CASE:
        for nu in NUS:
            for p in POS:
                c = cell(cl, z, nu, p); ent = []
                for k in RULES:
                    a = counts(c, k); fl = wilson(a[1], 1000)[1] > ALPHA
                    ent.append('%d/%d/%d%s' % (a[0], a[1], a[2], '$^\\dagger$' if fl else ''))
                s = counts(c, 'S'); cond = '%d/%d' % (s[1], s[0] + s[1]) if s[0] + s[1] else '0/0'
                g.write('%s $z{=}%g$ & %g & %s & $%s$ & %d & %s & %s\\\\\n' % (CN[cl], z, nu, PTEX[p], c['truth'], selC(c), ' & '.join(ent), cond))
    g.write('\\bottomrule\n\\end{longtable}\n')

# ---------------- Supplement: fitted-branch agreement per cell ----------------
with open(os.path.join(HERE, 'tableS_fit.tex'), 'w') as g:
    head = 'case & $\\nu$ & position & truth & fit outside & fit wrong & S$=$fit (fit wrong) & I$=$fit (fit wrong) & C$=$fit (fit wrong)\\\\\n'
    g.write('\\begin{longtable}{@{}lrllrrrrr@{}}\n\\caption{Agreement with the fitted optimizer, kept separate from branch declarations. ``S$=$fit\'\' counts draws in which S declares and the complex global fit lies in the declared arc; the parenthesis counts those in which the fitted branch is not the noiseless branch. The fits are numerical multistart optima.}\\label{tab:fit}\\\\\n\\toprule\n' + head + '\\midrule\n\\endfirsthead\n\\toprule\n' + head + '\\midrule\n\\endhead\n')
    for (cl, z) in CASE:
        for nu in NUS:
            for p in POS:
                c = cell(cl, z, nu, p); R = c['records']; tr = c['truth']
                ent = []
                for k in ('S', 'I', 'C'):
                    ag = [r for r in R if r['d' + k] != 'inc' and r['d' + k] == r['fit']]
                    ent.append('%d (%d)' % (len(ag), sum(r['fit'] != tr for r in ag)))
                g.write('%s $z{=}%g$ & %g & %s & $%s$ & %d & %d & %s\\\\\n' % (CN[cl], z, nu, PTEX[p], tr, sum(r['fit'] == 'out' for r in R),
                        sum(r['fit'] in ('A', 'B') and r['fit'] != tr for r in R), ' & '.join(ent)))
    g.write('\\bottomrule\n\\end{longtable}\n')

# ---------------- Supplement: scope table ----------------
SC = os.path.join(AS, 'scope', 'scope_cells.csv')
sc = [r for r in csv.DictReader(open(SC)) if r['study'].startswith('fresh')]
with open(os.path.join(HERE, 'tableS_scope.tex'), 'w') as g:
    g.write('\\begin{tabular}{@{}lccccc@{}}\n\\toprule\n & & & & certified & in certified\\\\\n case & cells & $|\\lambda/\\lambda_*-1|$ & $\\lambda_c/\\lambda_*-1$ & half-width$/\\lambda_*$ & interval\\\\\n\\midrule\n')
    for (cl, z) in CASE:
        rr = [r for r in sc if r['cluster'] == cl and float(r['z']) == z]
        dv = [abs(float(r['rel_dev'])) for r in rr]
        hw = '$5.8\\cdot10^{-5}$' if cl != 'pair' else 'not stated'
        st = 'none' if cl != 'pair' else 'undetermined'
        g.write('%s & %d & $%.1e$--$%.1e$ & $%+.2e$ & %s & %s\\\\\n' % (CTEX[(cl, z)], len(rr), min(dv), max(dv), float(rr[0]['lamc_rel_dev']), hw, st))
    g.write('\\bottomrule\n\\end{tabular}\n')

# summary numbers for the text
S_ = {k: (sum(counts(c, k)[0] for c in CELLS), sum(counts(c, k)[1] for c in CELLS)) for k in RULES}
print(json.dumps(S_), 'restored', sum(sum(r['dS'] == c['truth'] and r['dI'] == 'inc' for r in c['records']) for c in CELLS),
      'added wrong', sum(sum(r['dS'] not in ('inc', c['truth']) and r['dI'] == 'inc' for r in c['records']) for c in CELLS))
