"""Compact validation table: worst near-tie error (over tie-, tie+, between) per rule."""
import os, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, '..', 'results', 'validation_raw.json')))
from summarize import wilson
groups = {}
for c in D['configs']:
    groups.setdefault((c['cluster'], c['z'], c['nu']), {})[c['pos']] = c
lines = ['| cluster | z | ν | L | C | **I** | plug-in P | trust fit F | I certifies (far, 2 widths) | fit outside arcs (tie) |',
         '|---|---|---|---|---|---|---|---|---|---|']
for (n, z, nu), d in groups.items():
    def worst(rule):
        best = None
        for pos in ('tie-', 'tie+', 'between'):
            R = d[pos]['records']; tr = d[pos]['truth']
            if rule == 'F':
                k = sum(r['fit'] in ('A', 'B') and r['fit'] != tr for r in R)
            else:
                k = sum(r['d' + rule] not in ('inc', tr) for r in R)
            w = wilson(k, len(R))
            if best is None or w[0] > best[0]:
                best = w
        return '%.3f (%.3f–%.3f)' % best
    far = d['far+']['records']
    cert = sum(r['dI'] == r['fit'] and r['fit'] in ('A', 'B') for r in far) / len(far)
    out = sum(r['fit'] == 'out' for r in d['tie+']['records']) / len(d['tie+']['records'])
    lines.append(f'| {n} | {z} | {nu} | {worst("L")} | {worst("C")} | **{worst("I")}** | {worst("P")} | {worst("F")} | {cert:.2f} | {out:.2f} |')
txt = '\n'.join(lines)
open(os.path.join(HERE, '..', 'results', 'validation_table.md'), 'w').write(txt + '\n')
print(txt)
