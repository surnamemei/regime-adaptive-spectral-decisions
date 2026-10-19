"""Small targeted calibration study of the reliability rules (original model).

Clusters: equal pair (K=2), symmetric triple and asymmetric triple (K=3); N=21, b=10.
For each (cluster, z, noise level nu, lambda position) draw `ndraw` noise vectors and
record the decisions of rules L, C, I (ours), P (plug-in probit) and F (trust the fit).
Truth = pseudo-true label of the noiseless data (two Hermitian chart minima).
Usage: python validate.py [ndraw] [nworkers]
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import sys, json, math, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
from reliability import Setup, analyse, noiseless_label
from scipy.stats import t as tdist

OUT = os.path.join(HERE, '..', 'results')
N, B, ALPHA = 21, 10.0, 0.05
Q = float(tdist.ppf(1 - ALPHA, N))
S = Setup(N, B)
LAM = {1: 0.032990205563495685, 2: 2.5936704891874683e-4}
ARCS = {1: ((-8, -1), (10, 15)), 2: ((-15, -5), (12, 20))}
# reduced-problem constants (limit model, ../../results/weak_noise.json)
CONST = {1: dict(Phi1=1.114344120, sX=0.249435265, unit=0.0370563141),
         2: dict(Phi1=2.112466994e-3, sX=6.645863625e-4, unit=1.916883326e-4)}
CL = {'pair': ([1, 1], [-1, 1], 1, 2.0),
      'sym3': ([1, 1, 1], [-1, 0, 1], 2, 2 / 3),
      'asym3': ([1, 2, .5], [-1, .3, 1.5], 2, 1.9144115796098178)}
CASES = [('pair', 0.3), ('sym3', 0.3), ('asym3', 0.3), ('asym3', 0.15)]
NUS = [0.02, 0.1, 0.5, 1.0]
POS = ['tie-', 'tie+', 'between', 'far+']


def cz(w, o, z):
    return (np.array(w, float)[None, :] * np.exp(1j * z * np.outer(S.s, o))).sum(1)


def lam_c(name, z):
    w, o, r, h = CL[name]; A = h * z ** (2 * r); ls = LAM[r]
    D = lambda lam: noiseless_label(S, cz(w, o, z) - lam * A * np.exp(1j * B * S.s), r, ARCS[r])[1]
    lo, hi = ls * 0.95, ls * 1.05; dlo = D(lo)
    for _ in range(45):
        mid = (lo + hi) / 2; dm = D(mid)
        if (dm < 0) == (dlo < 0):
            lo, dlo = mid, dm
        else:
            hi = mid
    return (lo + hi) / 2


def job(args):
    name, z, nu, pos, lamc, seed, nd = args
    w, o, r, h = CL[name]; A = h * z ** (2 * r); c = CONST[r]
    rho = nu * c['unit']; sigma = rho * A; width = rho * c['sX'] / c['Phi1']; ls = LAM[r]
    lam = {'tie-': lamc - 0.1 * width, 'tie+': lamc + 0.1 * width,
           'between': 0.5 * (lamc + ls), 'far+': lamc + 2 * width}[pos]
    x = cz(w, o, z) - lam * A * np.exp(1j * B * S.s)
    truth, Dx = noiseless_label(S, x, r, ARCS[r])
    rng = np.random.default_rng(seed)
    recs = []
    for _ in range(nd):
        xi = (rng.normal(size=N) + 1j * rng.normal(size=N)) / math.sqrt(2)
        out = analyse(S, x + sigma * xi, r, ls, ARCS[r], Q, c['Phi1'], c['sX'], ALPHA)
        recs.append({k: out[k] for k in ('dL', 'dC', 'dI', 'dP', 'fit', 'tL', 'TC', 'A_hat', 'sigma_hat')})
    return dict(cluster=name, z=z, nu=nu, pos=pos, lam=lam, lam_c=lamc, lam_star=ls, A=A, sigma=sigma,
                rho=rho, width=width, truth=truth, D_noiseless=Dx, records=recs)


def main():
    nd = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    nw = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    chunk = 50
    lamcs = {(n, z): lam_c(n, z) for n, z in CASES}
    jobs = []
    seed = 20261020
    for n, z in CASES:
        for nu in NUS:
            for pos in POS:
                for k in range(nd // chunk):
                    seed += 1
                    jobs.append((n, z, nu, pos, lamcs[(n, z)], seed, chunk))
    t0 = time.time()
    with Pool(nw) as p:
        res = p.map(job, jobs, chunksize=1)
    merged = {}
    for rr in res:
        key = (rr['cluster'], rr['z'], rr['nu'], rr['pos'])
        if key not in merged:
            merged[key] = {k: v for k, v in rr.items() if k != 'records'}; merged[key]['records'] = []
        merged[key]['records'] += rr['records']
    out = dict(status='NUMERICAL_VALIDATION', N=N, b=B, alpha=ALPHA, t_quantile=Q, lambda_c=
               {f'{n}@{z}': v for (n, z), v in lamcs.items()}, seconds=time.time() - t0,
               configs=list(merged.values()))
    with open(os.path.join(OUT, 'validation_raw.json'), 'w') as f:
        json.dump(out, f)
    print('done', len(merged), time.time() - t0)


if __name__ == '__main__':
    main()
