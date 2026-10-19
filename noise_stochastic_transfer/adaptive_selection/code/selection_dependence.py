"""Supplementary analysis (not pre-specified; labelled exploratory in the report):
decomposition of S's wrong declarations by selected test, and dependence between
selection and each test's errors, in every cell with mixed selection."""
import os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, '..', 'results')
D = json.load(open(os.path.join(RES, 'raw.json')))
def wil(k, n, zq=1.959964):
    p = k / n; den = 1 + zq**2 / n; c = (p + zq**2/(2*n))/den; h = zq*math.sqrt(p*(1-p)/n + zq**2/(4*n*n))/den
    return p, max(0, c-h), min(1, c+h)
out = []
tot = dict(n=0, Sw=0, Cw=0, Lw=0, exceed=0, cells=0)
for c in D['configs']:
    R, tr = c['records'], c['truth']; n = len(R)
    selC = [r for r in R if r['sel'] == 'C']; selL = [r for r in R if r['sel'] == 'L']
    if not (20 <= len(selC) <= n - 20):
        continue
    w = lambda rs, k: sum(r['d' + k] not in ('inc', tr) for r in rs)
    Sw, Cw, Lw = w(R, 'S'), w(R, 'C'), w(R, 'L')
    tot['cells'] += 1; tot['n'] += n; tot['Sw'] += Sw; tot['Cw'] += Cw; tot['Lw'] += Lw; tot['exceed'] += Sw > max(Cw, Lw)
    # dependence of selection on the noise-scale factor V and on A_hat
    V = np.array([r['sigma_hat'] / c['sigma'] for r in R]); Ar = np.array([r['A_hat'] / c['A'] for r in R])
    s = np.array([r['sel'] == 'C' for r in R])
    out.append('%s z=%g nu=%g %s truth %s | selC %d | S wrong %d (via C %d, via L %d) | C wrong %d (|selC %d/%d, |selL %d/%d) | L wrong %d (|selL %d/%d, |selC %d/%d) | corr(sel,V) %.2f corr(sel,Ahat/A) %.2f' % (
        c['cluster'], c['z'], c['nu'], c['pos'], tr, len(selC), Sw, w(selC, 'C'), w(selL, 'L'),
        Cw, w(selC, 'C'), len(selC), w(selL, 'C'), len(selL), Lw, w(selL, 'L'), len(selL), w(selC, 'L'), len(selC),
        np.corrcoef(s, V)[0, 1], np.corrcoef(s, Ar)[0, 1]))
out.append('cells with mixed selection: %d, draws %d; pooled wrong S %d [%.4f, %.4f], C %d, L %d; cells with S wrong > max(C,L): %d' % (
    tot['cells'], tot['n'], tot['Sw'], *wil(tot['Sw'], tot['n'])[1:], tot['Cw'], tot['Lw'], tot['exceed']))
open(os.path.join(RES, 'selection_dependence.txt'), 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
