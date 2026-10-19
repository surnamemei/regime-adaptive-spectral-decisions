"""Exploratory summaries (not pre-specified): independent-selection benchmark in mixed cells,
spread of A_hat/A, and correct/wrong declarations by noise level and side (pooled, secondary)."""
import os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, '..', 'results')
D = json.load(open(os.path.join(RES, 'raw.json')))
L = []
bench_tot = obs_tot = 0
for c in D['configs']:
    R, tr = c['records'], c['truth']; n = len(R)
    selC = [r for r in R if r['sel'] == 'C']
    if 20 <= len(selC) <= n - 20:
        w = lambda k: sum(r['d' + k] not in ('inc', tr) for r in R)
        pc = len(selC) / n
        bench = pc * w('C') + (1 - pc) * w('L'); obs = w('S')
        bench_tot += bench; obs_tot += obs
L.append('mixed-selection cells: observed S wrong %d vs independent-selection benchmark %.1f (ratio %.3f)' % (obs_tot, bench_tot, obs_tot / bench_tot))
for nu in sorted({c['nu'] for c in D['configs']}):
    cs = [c for c in D['configs'] if c['nu'] == nu]
    ar = np.concatenate([[r['A_hat'] / c['A'] for r in c['records']] for c in cs])
    selC = sum(r['sel'] == 'C' for c in cs for r in c['records']); n = sum(len(c['records']) for c in cs)
    line = 'nu=%g: sd(A_hat/A) %.3f, P(A_hat<=0) %.4f, selected C %d/%d |' % (nu, ar.std(), np.mean(ar <= 0), selC, n)
    for side in ('A', 'B'):
        cc = [c for c in cs if c['truth'] == side]
        for k in ('S', 'I', 'C', 'L', 'C2'):
            cor = sum(r['d' + k] == c['truth'] for c in cc for r in c['records'])
            wr = sum(r['d' + k] not in ('inc', c['truth']) for c in cc for r in c['records'])
            line += ' %s:%s %d/%d' % (side, k, cor, wr)
        line += ' |'
    L.append(line)
open(os.path.join(RES, 'gains_summary.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
