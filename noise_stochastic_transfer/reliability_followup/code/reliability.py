"""Reliability rules for the deficient-order fit (numerical implementation).

Observation: y = c_z - eps e_b + sigma xi (xi ~ CN(0, I)), centred grid s = t/N.
Known: N, K = r+1, b, the reduced threshold lambda_* and the wells U_A, U_B.
Unknown: the cluster, z, eps, sigma, hence A = h_r z^{2r}, lambda = eps/A, rho = sigma/A.

Statistics (all computed from y only):
  sigma_hat^2 = 2 ||iota y_A||^2 / N            (anti-Hermitian channel; chi^2_N exact)
  t_L = (eps_hat - lambda_* A_hat) / se          (Gauss-defect linear contrast)
  T_C = (J_A - J_B) / (sqrt2 sigma_hat ||r_A - r_B||)   (Hermitian chart costs)
Rules: L, C, and I = declare a label only if L and C declare the same label.
Reliability output: 'reliable' iff the rule declares X and the complex global fit lies in U_X.
Baselines: trust-the-fit; plug-in probit.
NUMERICAL CODE; the calibration claim is the theorem in ../RELIABILITY_RESULT.md.
"""
import math
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm, t as tdist

N_DEFAULT = 21


def sig(n):
    return 1 if n % 4 in (0, 1) else -1


class Setup:
    def __init__(self, N=21, b=10.0):
        self.N = N
        self.T = (N - 1) // 2
        self.s = np.arange(-self.T, self.T + 1) / N
        self.b = b

    # ---- Hermitian / anti-Hermitian channels -------------------------------
    def split(self, y):
        yH = (y + np.conj(y[::-1])) / 2
        yA = (y - np.conj(y[::-1])) / 2j
        return yH.real + yH.imag, yA.real + yA.imag

    def psi(self, u):
        a = np.outer(self.s, np.atleast_1d(u))
        return np.cos(a) + np.sin(a)

    def dpsi(self, u):
        a = np.outer(self.s, np.atleast_1d(u))
        return self.s[:, None] * (np.cos(a) - np.sin(a))


# ---------------------------------------------------------------------------
# Real (Hermitian-channel) variable projection
class HermLS:
    def __init__(self, S, h):
        self.S, self.h = S, h

    def resid(self, u):
        Q, _ = np.linalg.qr(self.S.psi(u))
        return self.h - Q @ (Q.T @ self.h)

    def jac(self, u):
        E = self.S.psi(u)
        Q, R = np.linalg.qr(E)
        a = np.linalg.solve(R, Q.T @ self.h)
        r = self.h - Q @ (Q.T @ self.h)
        D = self.S.dpsi(u)
        RinvT = np.linalg.inv(R).T
        J = np.empty((len(self.h), len(u)))
        for k in range(len(u)):
            t1 = D[:, k] * a[k]
            t1 = t1 - Q @ (Q.T @ t1)
            ek = np.zeros(len(u)); ek[k] = D[:, k] @ r
            J[:, k] = -(t1 + Q @ (RinvT @ ek))
        return J

    def fit(self, u0, free_idx, lo=None, hi=None):
        u0 = np.asarray(u0, float)
        fi = list(free_idx)

        def full(p):
            u = u0.copy(); u[fi] = p; return u
        lb = np.full(len(fi), -np.inf); ub = np.full(len(fi), np.inf)
        if lo is not None:
            # last free index is the free frequency if bounded
            lb[-1], ub[-1] = lo, hi
        p0 = np.clip(u0[fi], lb + 1e-12, ub - 1e-12)
        res = least_squares(lambda p: self.resid(full(p)), p0,
                            jac=lambda p: self.jac(full(p))[:, fi],
                            method='trf' if lo is not None else 'lm', bounds=(lb, ub) if lo is not None else (-np.inf, np.inf),
                            xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
        u = full(res.x)
        r = self.resid(u)
        return u, float(r @ r), r


def local_init(S, h, r):
    """Data-driven start for the r local nodes: best r-term Hermitian fit near 0."""
    H = HermLS(S, h)
    best = None
    for c in [0.01, 0.03, 0.1, 0.3, 1.0]:
        u0 = np.array([0.0]) if r == 1 else np.linspace(-c, c, r)
        if r == 1:
            u0 = np.array([c * 0.1])
        u, cost, _ = H.fit(u0, range(r))
        if best is None or cost < best[1]:
            best = (u, cost)
    return np.sort(best[0])


def chart_min(S, h, r, uloc, arc, step=0.25):
    """Hermitian-channel minimum over K-term fits with one root in the arc."""
    H = HermLS(S, h)
    lo, hi = arc
    grid = np.arange(lo, hi + 1e-9, step)
    loc = uloc.copy()
    prof = []
    for v in grid:
        u, c, _ = H.fit(np.r_[loc, v], range(r))
        loc = u[:r]
        prof.append((c, u.copy()))
    order = np.argsort([p[0] for p in prof])
    best = None
    for k in order[:3]:
        u0 = prof[k][1]
        u, c, res = H.fit(u0, range(r + 1), lo=lo, hi=hi)
        if best is None or c < best[0]:
            best = (c, u, res)
    return best


# ---------------------------------------------------------------------------
# Complex global fit: label of the global optimum
def complex_global_label(S, y, r, uloc, arcs, step=0.5):
    from original_fit import Problem  # noqa (path set by caller)
    N = S.N
    K = r + 1
    pb = Problem.__new__(Problem)
    pb.N, pb.T, pb.s, pb.K, pb.y = N, S.T, S.s, K, y
    grid = np.arange(-math.pi * N, math.pi * N, step)
    loc = uloc.copy()
    prof = []
    for v in grid:
        u, c = pb.refine(np.r_[loc, v], fixed={K - 1: v})
        if np.max(np.abs(u[:r])) > 3.0:
            u, c = pb.refine(np.r_[uloc, v], fixed={K - 1: v})
        loc = u[:r]
        prof.append((c, u.copy()))
    idx = [i for i in range(len(grid)) if prof[i][0] <= prof[i - 1][0] and prof[i][0] <= prof[(i + 1) % len(grid)][0]]
    idx = sorted(idx, key=lambda i: prof[i][0])[:4]
    cands = []
    for i in idx:
        u, c = pb.refine(prof[i][1])
        cands.append((c, u))
    for sc in [0.5, 1.0, 2.0]:
        u, c = pb.refine(np.r_[uloc, sc * (uloc.max() + 0.1)])
        cands.append((c, u))
    c, u = min(cands, key=lambda x: x[0])
    P = 2 * math.pi * N
    uw = (u + P / 2) % P - P / 2
    v = uw[np.argmax(np.abs(uw))]
    if abs(v) < 1.0:
        return 'out', float(v)
    for lab, (lo, hi) in zip(('A', 'B'), arcs):
        if lo <= v <= hi:
            return lab, float(v)
    return 'out', float(v)


# ---------------------------------------------------------------------------
# Linear Gauss-defect statistic
def gauss_defect(m, r):
    H = np.array([[m[i + j] for j in range(r)] for i in range(r)])
    bv = np.array([m[r + i] for i in range(r)])
    return m[2 * r] - bv @ np.linalg.solve(H, bv)


def linear_stat(S, h, r, lam_star, sigma_hat):
    X = np.column_stack([S.s ** n for n in range(2 * r + 1)] + [S.psi(S.b)[:, 0]])
    XtX_inv = np.linalg.inv(X.T @ X)
    theta = XtX_inv @ (X.T @ h)

    def f(th):
        m = np.array([math.factorial(n) * sig(n) * th[n] for n in range(2 * r + 1)])
        A = gauss_defect(m, r)
        eps = -th[-1]
        return eps - lam_star * A, A, eps

    fval, A_hat, eps_hat = f(theta)
    g = np.empty(len(theta))
    for k in range(len(theta)):
        dk = 1e-6 * max(abs(theta[k]), 1e-8)
        tp = theta.copy(); tp[k] += dk
        tm = theta.copy(); tm[k] -= dk
        g[k] = (f(tp)[0] - f(tm)[0]) / (2 * dk)
    se = sigma_hat * math.sqrt(g @ XtX_inv @ g / 2)
    return fval / se, A_hat, eps_hat


# ---------------------------------------------------------------------------
def decide(stat, q):
    if stat < -q:
        return 'A'
    if stat > q:
        return 'B'
    return 'inc'


def analyse(S, y, r, lam_star, arcs, q, Phi1, sX, alpha):
    """All statistics and decisions for one data vector."""
    h, a = S.split(y)
    N = S.N
    sigma_hat = math.sqrt(2 * (a @ a) / N)
    uloc = local_init(S, h, r)
    JA, uA, rA = chart_min(S, h, r, uloc, arcs[0])
    JB, uB, rB = chart_min(S, h, r, uloc, arcs[1])
    TC = (JA - JB) / (math.sqrt(2) * sigma_hat * np.linalg.norm(rA - rB))
    tL, A_hat, eps_hat = linear_stat(S, h, r, lam_star, sigma_hat)
    dL, dC = decide(tL, q), decide(TC, q)
    dI = dL if dL == dC else 'inc'
    fit_label, vfit = complex_global_label(S, y, r, uloc, arcs)
    # plug-in probit baseline (predictive probability that the fit selects A)
    if A_hat > 0:
        lam_hat, rho_hat = eps_hat / A_hat, sigma_hat / A_hat
        pA = float(norm.cdf(-Phi1 * (lam_hat - lam_star) / (rho_hat * sX)))
        dP = 'A' if pA >= 1 - alpha else ('B' if pA <= alpha else 'inc')
    else:
        pA, dP = float('nan'), 'inc'
    return dict(sigma_hat=sigma_hat, tL=float(tL), TC=float(TC), A_hat=float(A_hat), eps_hat=float(eps_hat),
                dL=dL, dC=dC, dI=dI, dP=dP, pA=pA, fit=fit_label, vfit=vfit,
                JA=JA, JB=JB, vA=float(uA[-1]), vB=float(uB[-1]))


def noiseless_label(S, x_complex, r, arcs):
    h, _ = S.split(x_complex)
    uloc = local_init(S, h, r)
    JA = chart_min(S, h, r, uloc, arcs[0])[0]
    JB = chart_min(S, h, r, uloc, arcs[1])[0]
    return ('A' if JA < JB else 'B'), JA - JB
