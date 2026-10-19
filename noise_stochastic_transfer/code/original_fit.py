"""Direct global fits of the ORIGINAL model (numerical, not certified).

Data  y = c_z - eps e_b + sigma xi  on s = t/N, fitted in l2 by K exponentials
with real frequencies and complex amplitudes (variable projection).
Global search does not presuppose the allocation found by the theory:
  (i)   free-frequency profile: K-1 local nodes optimized for each v on a grid,
        then full K-dimensional refinement of the best profile minima;
  (ii)  all-local starts (all K nodes near 0, several scales);
  (iii) uniformly random starts on the frequency torus.
Final candidates can be re-evaluated and Gauss-Newton refined in mpmath.
"""
import math
import numpy as np
from scipy.optimize import least_squares


class Problem:
    def __init__(self, N, weights, offsets, K, z, eps, sigma, xi, b=10.0):
        self.N = N
        self.T = (N - 1) // 2
        self.s = np.arange(-self.T, self.T + 1) / N
        self.K = K
        self.z = z
        self.weights = np.asarray(weights, float)
        self.offsets = np.asarray(offsets, float)
        cz = (self.weights[None, :] * np.exp(1j * z * np.outer(self.s, self.offsets))).sum(1)
        self.y = cz - eps * np.exp(1j * b * self.s) + sigma * xi
        self.ynorm = np.linalg.norm(self.y)

    def E(self, u):
        return np.exp(1j * np.outer(self.s, u))

    def resid(self, u):
        Q, _ = np.linalg.qr(self.E(u))
        return self.y - Q @ (Q.conj().T @ self.y)

    def cost(self, u):
        r = self.resid(u)
        return float(np.vdot(r, r).real)

    def amplitudes(self, u):
        return np.linalg.lstsq(self.E(u), self.y, rcond=None)[0]

    def _rvec(self, u):
        r = self.resid(u)
        return np.r_[r.real, r.imag]

    def _jac(self, u):
        """Exact Golub-Pereyra Jacobian of the VarPro residual (I-P)y in u."""
        E = self.E(u)
        Q, R = np.linalg.qr(E)
        a = np.linalg.solve(R, Q.conj().T @ self.y)
        r = self.y - Q @ (Q.conj().T @ self.y)
        J = np.empty((len(self.s), len(u)), complex)
        Rinvh = np.linalg.inv(R).conj().T
        for k in range(len(u)):
            Dk = 1j * self.s * E[:, k]
            t1 = Dk * a[k]
            t1 = t1 - Q @ (Q.conj().T @ t1)
            ek = np.zeros(len(u), complex)
            ek[k] = np.vdot(Dk, r)
            t2 = Q @ (Rinvh @ ek)
            J[:, k] = -(t1 + t2)
        return np.r_[J.real, J.imag]

    def refine(self, u0, fixed=None, scale=None):
        """Local least-squares refinement. fixed: dict index->value (held).
        scale: per-parameter scale (local nodes are optimized in units of z)."""
        u0 = np.asarray(u0, float)
        free = [k for k in range(len(u0)) if not fixed or k not in fixed]
        sc = np.ones(len(u0)) if scale is None else np.asarray(scale, float)

        def full(p):
            u = u0.copy()
            u[free] = p * sc[free]
            if fixed:
                for k, val in fixed.items():
                    u[k] = val
            return u

        p0 = u0[free] / sc[free]
        jac = lambda p: self._jac(full(p))[:, free] * sc[free][None, :]
        res = least_squares(lambda p: self._rvec(full(p)), p0, jac=jac, method='lm',
                            xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        u = full(res.x)
        return u, self.cost(u)


def wrap(v, N):
    P = 2 * math.pi * N
    return (v + P / 2) % P - P / 2


def global_fit(pb, gauss_nodes, rng, vstep=0.25, n_random=60, n_top=6, verbose=False):
    N, K, z = pb.N, pb.K, pb.z
    r = K - 1
    scale_loc = np.r_[np.full(r, z), 1.0]
    cands = []
    # (i) free-frequency profile with optimized local nodes (warm-started)
    vgrid = np.arange(-math.pi * N, math.pi * N, vstep)
    prof = np.empty(len(vgrid))
    loc = z * np.asarray(gauss_nodes, float)
    locs = []
    for i, v in enumerate(vgrid):
        u0 = np.r_[loc, v]
        u, c = pb.refine(u0, fixed={K - 1: v}, scale=scale_loc)
        # guard against drift of a local node far away: reset if needed
        if np.max(np.abs(u[:r])) > 20 * z * (1 + np.max(np.abs(gauss_nodes))):
            u, c = pb.refine(np.r_[z * np.asarray(gauss_nodes), v], fixed={K - 1: v}, scale=scale_loc)
        prof[i] = c
        loc = u[:r]
        locs.append(u[:r].copy())
    # local minima of the profile (circular)
    idx = [i for i in range(len(vgrid)) if prof[i] <= prof[i - 1] and prof[i] <= prof[(i + 1) % len(vgrid)]]
    idx = sorted(idx, key=lambda i: prof[i])[:n_top]
    for i in idx:
        u, c = pb.refine(np.r_[locs[i], vgrid[i]], scale=scale_loc)
        cands.append(('profile', u, c))
    # (ii) all-local starts
    base = z * np.asarray(gauss_nodes, float)
    for sc in [z, 3 * z, math.sqrt(z), 0.3, 1.0, 2.0]:
        for zeta in [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]:
            u, c = pb.refine(np.r_[base, zeta * sc], scale=np.r_[np.full(r, z), sc])
            cands.append(('alllocal', u, c))
    # (iii) random starts
    for _ in range(n_random):
        u0 = rng.uniform(-math.pi * N, math.pi * N, size=K)
        u, c = pb.refine(u0)
        cands.append(('random', u, c))
    cands.sort(key=lambda x: x[2])
    best = cands[0]
    return best, cands, (vgrid, prof)


def mp_refine(N, weights, offsets, z, eps, sigma, xi, u, b=10.0, dps=40, iters=80):
    """Gauss-Newton refinement of the frequencies in mpmath with exact
    (Golub-Pereyra) VarPro Jacobian. The data are regenerated in mpmath from
    the same float noise vector xi, so no float64 cancellation enters.
    Returns (u, cost, grad_norm) with the cost as an mpf."""
    import mpmath as mp
    mp.mp.dps = dps
    T = (N - 1) // 2
    s = [mp.mpf(t) / N for t in range(-T, T + 1)]
    z, eps, sigma, b = mp.mpf(z), mp.mpf(eps), mp.mpf(sigma), mp.mpf(b)
    y = mp.matrix([sum(mp.mpf(w) * mp.expj(z * mp.mpf(o) * sj) for w, o in zip(weights, offsets))
                   - eps * mp.expj(b * sj) + sigma * mp.mpc(complex(x)) for sj, x in zip(s, xi)])
    u = [mp.mpf(float(x)) for x in u]
    K = len(u)

    def parts(u):
        E = mp.matrix([[mp.expj(uk * sj) for uk in u] for sj in s])
        G = E.H * E
        a = mp.lu_solve(G, E.H * y)
        r = y - E * a
        return E, G, a, r

    def jac(E, G, a, r, u):
        cols = []
        Ginv = mp.inverse(G)
        for k in range(K):
            Dk = mp.matrix([1j * sj * E[i, k] for i, sj in enumerate(s)])
            t1 = Dk * a[k]
            t1 = t1 - E * (Ginv * (E.H * t1))
            ek = mp.matrix(K, 1)
            ek[k] = (Dk.H * r)[0]
            t2 = E * (Ginv * ek)
            cols.append(-(t1 + t2))
        return cols

    def system(u):
        E, G, a, r = parts(u)
        cols = jac(E, G, a, r, u)
        n = len(s)
        J = mp.matrix(2 * n, K)
        for k in range(K):
            for i in range(n):
                J[i, k] = mp.re(cols[k][i])
                J[n + i, k] = mp.im(cols[k][i])
        rv = mp.matrix([mp.re(x) for x in r] + [mp.im(x) for x in r])
        return J, rv, sum(x ** 2 for x in rv)

    mu = mp.mpf('1e-6')
    J, rv, c0 = system(u)
    for _ in range(iters):
        JtJ, g = J.T * J, J.T * rv
        diag = [JtJ[k, k] for k in range(K)]
        while True:
            A = JtJ.copy()
            for k in range(K):
                A[k, k] += mu * diag[k]
            step = mp.lu_solve(A, -g)
            un = [u[k] + step[k] for k in range(K)]
            Jn, rvn, cn = system(un)
            if cn <= c0:
                u, J, rv, c0 = un, Jn, rvn, cn
                mu = max(mu / 10, mp.mpf('1e-30'))
                break
            mu *= 10
            if mu > 1e20:
                break
        gnorm = max(abs(x) for x in (J.T * rv))
        if gnorm < mp.mpf(10) ** (-(dps - 8)) * max(c0, mp.mpf(10) ** (-dps)):
            break
    E, G, a, r = parts(u)
    cols = jac(E, G, a, r, u)
    grad = [2 * sum(mp.re(mp.conj(r[i]) * cols[k][i]) for i in range(len(s))) for k in range(K)]
    cost = sum(abs(x) ** 2 for x in r)
    return [float(x) for x in u], cost, float(max(abs(g) for g in grad)), [complex(x) for x in a]
