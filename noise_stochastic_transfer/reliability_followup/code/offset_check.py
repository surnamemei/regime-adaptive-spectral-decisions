"""Noiseless offset of the linear statistic: (f_hat(x) - (eps - lam_c A)) in units of
the L noise width at each noise level (explains where rule L alone is miscalibrated)."""
import sys, os, math, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from reliability import Setup, linear_stat
import validate as V
S = V.S
out = []
D = json.load(open(os.path.join(HERE, '..', 'results', 'validation_raw.json')))
for n, z in V.CASES:
    w, o, r, h = V.CL[n]; A = h * z ** (2 * r); ls = V.LAM[r]; lc = D['lambda_c'][f'{n}@{z}']
    lam = lc
    x = V.cz(w, o, z) - lam * A * np.exp(1j * V.B * S.s)
    hH, _ = S.split(x)
    # unit sigma_hat to get f_hat and se per unit sigma
    t1, Ahat, ehat = linear_stat(S, hH, r, ls, 1.0)
    fhat = ehat - ls * Ahat
    se1 = fhat / t1                       # se for sigma=1
    for nu in V.NUS:
        sigma = nu * V.CONST[r]['unit'] * A
        out.append(dict(cluster=n, z=z, nu=nu, offset_in_L_sd=fhat / (se1 * sigma),
                        Ahat_rel_bias=Ahat / A - 1, eps_hat_rel_bias=ehat / (lam * A) - 1,
                        lamc_shift_rel=lc / ls - 1))
        print(out[-1])
json.dump(out, open(os.path.join(HERE, '..', 'results', 'L_offset.json'), 'w'), indent=1)
