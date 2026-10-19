"""Summarize results/validation_raw.json: error rates with Wilson 95% intervals."""
import os, json, math, csv
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, '..', 'results', 'validation_raw.json')))


def wilson(k, n, zq=1.959964):
    if n == 0:
        return (float('nan'),) * 3
    p = k / n
    den = 1 + zq ** 2 / n
    c = (p + zq ** 2 / (2 * n)) / den
    h = zq * math.sqrt(p * (1 - p) / n + zq ** 2 / (4 * n ** 2)) / den
    return p, max(0.0, c - h), min(1.0, c + h)


def main():
    rows = []
    for cfg in D['configs']:
        R = cfg['records']; n = len(R); tr = cfg['truth']
        row = dict(cluster=cfg['cluster'], z=cfg['z'], nu=cfg['nu'], pos=cfg['pos'], truth=tr, n=n,
                   shift_over_width=(cfg['lam'] - cfg['lam_star']) / cfg['width'],
                   lamc_shift_over_width=(cfg['lam_c'] - cfg['lam_star']) / cfg['width'])
        row['P_fit_out'] = sum(r['fit'] == 'out' for r in R) / n
        for rule in ('L', 'C', 'I', 'P'):
            k = sum(r['d' + rule] not in ('inc', tr) for r in R)
            row[f'wrong_{rule}'] = wilson(k, n)
            kr = sum(r['d' + rule] == r['fit'] and r['fit'] in ('A', 'B') and r['fit'] != tr for r in R)
            row[f'falserel_{rule}'] = wilson(kr, n)
            row[f'reliable_{rule}'] = sum(r['d' + rule] == r['fit'] and r['fit'] in ('A', 'B') for r in R) / n
        kf = sum(r['fit'] in ('A', 'B') and r['fit'] != tr for r in R)
        row['falserel_F'] = wilson(kf, n)
        row['reliable_F'] = sum(r['fit'] in ('A', 'B') for r in R) / n
        row['P_Ahat_nonpos'] = sum(r['A_hat'] <= 0 for r in R) / n
        rows.append(row)

    fmt = lambda t: '%.3f [%.3f,%.3f]' % t
    with open(os.path.join(HERE, '..', 'results', 'validation_summary.csv'), 'w', newline='') as f:
        keys = list(rows[0].keys())
        w = csv.writer(f); w.writerow(keys)
        for r in rows:
            w.writerow([fmt(r[k]) if isinstance(r[k], tuple) else r[k] for k in keys])
    print('cluster z nu pos truth | wrong L | wrong C | wrong I | wrong P || falseRel I | falseRel P | falseRel F || rel I | rel F | out')
    for r in rows:
        print('%-5s %.2f %.2f %-7s %s | %s | %s | %s | %s || %s | %s | %s || %.3f | %.3f | %.3f' % (
            r['cluster'], r['z'], r['nu'], r['pos'], r['truth'], fmt(r['wrong_L']), fmt(r['wrong_C']), fmt(r['wrong_I']),
            fmt(r['wrong_P']), fmt(r['falserel_I']), fmt(r['falserel_P']), fmt(r['falserel_F']), r['reliable_I'],
            r['reliable_F'], r['P_fit_out']))


if __name__ == '__main__':
    main()
