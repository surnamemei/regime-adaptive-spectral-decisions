"""Fresh evaluation study of the regime-selection rule (PROTOCOL.md).

Reuses the existing implementation in reliability_followup/code unchanged and adds:
  - the selector s(y) of THEOREM.md section 1 (gamma = 1/(4r), c_nu = noise unit);
  - the fixed stricter-threshold cost test (level alpha/2);
  - per-draw records of all decisions, the selector inputs and the fitted branch.
Usage: python run_study.py [ndraw_per_cell] [nworkers] [--smoke]
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import sys, json, math, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
OLD = os.path.join(HERE, '..', '..', 'reliability_followup', 'code')
sys.path.insert(0, OLD); sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
from reliability import Setup, analyse, noiseless_label, gauss_defect, sig  # noqa: E402
from scipy.stats import t as tdist  # noqa: E402

N, B, ALPHA = 21, 10.0, 0.05
Q = float(tdist.ppf(1 - ALPHA, N))
Q2 = float(tdist.ppf(1 - ALPHA / 2, N))
S = Setup(N, B)
LAM = {1: 0.032990205563495685, 2: 2.5936704891874683e-4}
ARCS = {1: ((-8, -1), (10, 15)), 2: ((-15, -5), (12, 20))}
CONST = {1: dict(Phi1=1.114344120, sX=0.249435265, unit=0.0370563141),
         2: dict(Phi1=2.112466994e-3, sX=6.645863625e-4, unit=1.916883326e-4)}
CL = {'pair': ([1, 1], [-1, 1], 1, 2.0),
      'sym3': ([1, 1, 1], [-1, 0, 1], 2, 2 / 3),
      'asym3': ([1, 2, .5], [-1, .3, 1.5], 2, 1.9144115796098178)}
CASES = [('pair', 0.3), ('sym3', 0.3), ('asym3', 0.3), ('asym3', 0.15)]
NUS = [0.005, 0.02, 0.1, 0.35, 0.5, 0.7, 1.0]
POS = ['m2', 'm0.5', 'm0.1', 'between', 'p0.1', 'p0.5', 'p2']
XPOS = {'m2': -2.0, 'm0.5': -0.5, 'm0.1': -0.1, 'p0.1': 0.1, 'p0.5': 0.5, 'p2': 2.0}
GAMMA = {1: 1 / 4, 2: 1 / 8}          # gamma = 1/(4r)
MASTER_SEED = 2026101902
# lambda_c(z) from the earlier noiseless bisection (deterministic; reliability_followup/results)
LAMC = json.load(open(os.path.join(HERE, '..', '..', 'reliability_followup', 'results', 'validation_raw.json')))['lambda_c']


def cz(w, o, z):
    return (np.array(w, float)[None, :] * np.exp(1j * z * np.outer(S.s, o))).sum(1)


def selector(h, r, sigma_hat):
    """THEOREM.md section 1: returns (A_hat, m0_hat, nu_hat, S_hat, sel)."""
    X = np.column_stack([S.s ** n for n in range(2 * r + 1)] + [S.psi(S.b)[:, 0]])
    theta = np.linalg.solve(X.T @ X, X.T @ h)
    m = np.array([math.factorial(n) * sig(n) * theta[n] for n in range(2 * r + 1)])
    try:
        A_hat = float(gauss_defect(m, r))
    except np.linalg.LinAlgError:
        A_hat = 0.0
    m0 = float(m[0])
    unit = CONST[r]['unit']
    if A_hat > 0 and m0 > 0:
        nu_hat = sigma_hat / (unit * A_hat); S_hat = A_hat / m0
        sel = 'C' if nu_hat <= S_hat ** GAMMA[r] else 'L'
    else:
        nu_hat, S_hat, sel = float('nan'), float('nan'), 'L'
    return A_hat, m0, nu_hat, S_hat, sel


def decide(stat, q):
    return 'A' if stat < -q else ('B' if stat > q else 'inc')


def job(args):
    name, z, nu, pos, seed_words, nd = args
    w, o, r, h = CL[name]; A = h * z ** (2 * r); c = CONST[r]
    rho = nu * c['unit']; sigma = rho * A; width = rho * c['sX'] / c['Phi1']; ls = LAM[r]
    lamc = LAMC[f'{name}@{z}']
    lam = 0.5 * (lamc + ls) if pos == 'between' else lamc + XPOS[pos] * width
    x = cz(w, o, z) - lam * A * np.exp(1j * B * S.s)
    truth, Dx = noiseless_label(S, x, r, ARCS[r])
    rng = np.random.default_rng(np.random.SeedSequence(seed_words))
    recs = []
    for _ in range(nd):
        xi = (rng.normal(size=N) + 1j * rng.normal(size=N)) / math.sqrt(2)
        y = x + sigma * xi
        out = analyse(S, y, r, ls, ARCS[r], Q, c['Phi1'], c['sX'], ALPHA)
        hH, _ = S.split(y)
        A_hat, m0, nu_hat, S_hat, sel = selector(hH, r, out['sigma_hat'])
        dS = out['dC'] if sel == 'C' else out['dL']
        recs.append(dict(dL=out['dL'], dC=out['dC'], dI=out['dI'], dS=dS, dC2=decide(out['TC'], Q2),
                         sel=sel, fit=out['fit'], tL=out['tL'], TC=out['TC'], A_hat=A_hat, A_hat_L=out['A_hat'],
                         m0_hat=m0, sigma_hat=out['sigma_hat'], nu_hat=nu_hat, S_hat=S_hat))
    return dict(cluster=name, z=z, nu=nu, pos=pos, lam=lam, lam_c=lamc, lam_star=ls, A=A, sigma=sigma, rho=rho,
                width=width, truth=truth, D_noiseless=Dx, mu0=float(sum(w)), records=recs)


def main():
    smoke = '--smoke' in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    nd = int(args[0]) if args else 1000
    nw = int(args[1]) if len(args) > 1 else 30
    chunk = 50 if not smoke else nd
    jobs = []
    cells = [(n, z, nu, pos) for n, z in CASES for nu in NUS for pos in POS]
    if smoke:
        cells = [cells[0], cells[-1]]
    for ci, (n, z, nu, pos) in enumerate(cells):
        for k in range(max(1, nd // chunk)):
            jobs.append((n, z, nu, pos, [MASTER_SEED if not smoke else 1, ci, k], chunk))
    t0 = time.time()
    with Pool(nw) as p:
        res = p.map(job, jobs, chunksize=1)
    merged = {}
    for rr in res:
        key = (rr['cluster'], rr['z'], rr['nu'], rr['pos'])
        if key not in merged:
            merged[key] = {k: v for k, v in rr.items() if k != 'records'}; merged[key]['records'] = []
        merged[key]['records'] += rr['records']
    out = dict(status='SMOKE_TEST_DISCARD' if smoke else 'PROTOCOL_EVALUATION', N=N, b=B, alpha=ALPHA, q=Q, q_half=Q2,
               gamma=GAMMA, master_seed=MASTER_SEED, seconds=time.time() - t0, configs=list(merged.values()))
    fn = 'smoke_DISCARD.json' if smoke else 'raw.json'
    with open(os.path.join(HERE, '..', 'results', fn), 'w') as f:
        json.dump(out, f)
    print('done', len(merged), 'cells', time.time() - t0, 's')


if __name__ == '__main__':
    main()
