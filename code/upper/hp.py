"""High-precision position-side Weil form on the first K Legendre modes of one parity (ball arithmetic, python-flint).
Q_ab = Int_0^2 G_ab(r) dnu_L(r) + (d_L + c_L) delta_ab, with r = x/L and G_ab(r) = g_ab(Lr) the overlap of e_a with e_b shifted by x
(independent of L, an exact polynomial of degree <= D). So the integral is exact as sum_i w_i(L) G(r_i) over D+1 Chebyshev nodes,
with w_i built from the Chebyshev moments of the measure nu_L = pole/cosh weight + archimedean weight - prime spikes.
Same formula as decomp.pieces (float64), which matches the certified matrices."""
import time, numpy as np
from flint import arb, arb_mat, ctx

def gl(n):  # Gauss-Legendre nodes and weights on [-1, 1]
    return [arb.legendre_p_root(n, k, weight=True) for k in range(n)]

def legendre_rows(ts, kmax):  # P_0..P_kmax at each t (object arrays)
    P = [np.ones(len(ts), dtype=object) * arb(1), np.array(ts, dtype=object)]
    for k in range(1, kmax):
        P.append(((2 * k + 1) * ts * P[k] - k * P[k - 1]) / (k + 1))
    return P

class Overlap:
    """G(r_i) at the D+1 Chebyshev nodes r_i = 1 + u_i, for modes k = 2m + par, m < K."""
    def __init__(self, K, par):
        self.K, self.par = K, par
        kmax = 2 * K - 2 + par; ks = range(par, kmax + 1, 2)
        self.D = 2 * kmax + 1
        n = self.D + 1
        self.u = [arb.cos_pi(arb(2 * i + 1) / (2 * n)) for i in range(n)]
        g = gl(kmax + 2)  # exact for polynomial degree <= 2 kmax + 3 in t
        s = np.array([x for x, _ in g], dtype=object); w = np.array([y for _, y in g], dtype=object)
        nrm = [((2 * k + 1) * arb(1) / 2).sqrt() for k in ks]
        self.G = []
        for u in self.u:
            r = u + 1; h = (2 - r) / 2; t = -1 + h * (s + 1)
            P0, P1 = legendre_rows(t, kmax), legendre_rows(t + r, kmax)
            E0 = arb_mat([list(nrm[i] * P0[k] * w * h) for i, k in enumerate(ks)])
            E1 = arb_mat([list(nrm[i] * P1[k]) for i, k in enumerate(ks)])
            M = E0 * E1.transpose()
            self.G.append((M + M.transpose()) / 2)

def primes_upto(X):
    s = np.ones(X + 1, bool); s[:2] = False
    for i in range(2, int(X ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def cheb_moments(xs, c, D):  # sum_n c_n T_j(x_n), j = 0..D (object arrays of arb)
    T0 = np.ones(len(xs), dtype=object) * arb(1); T1 = xs.copy()
    out = [c.sum(), (c * T1).sum()]
    for _ in range(2, D + 1):
        T0, T1 = T1, 2 * xs * T1 - T0
        out.append((c * T1).sum())
    return out

def weil_matrix(O, L, Mq=None, parts=False):
    """Q (K x K arb_mat) at window half-width L (an arb or string)."""
    L = arb(L); D = O.D; Mq = Mq or 2 * D + 600  # arch weight has poles at r = i pi m / L; quadrature error is not in the radii (review)
    q = gl(Mq); r = np.array([x + 1 for x, _ in q], dtype=object); w = np.array([y for _, y in q], dtype=object)
    e = lambda z: np.array([arb.exp(v) for v in z], dtype=object)
    cosh_w = w * 4 * L * (e(L * r / 2) + e(-L * r / 2)) / 2
    k_w = w * 2 * L * e(-L * r / 2) / (1 - e(-2 * L * r))
    d = (w * 2 * L * (e(-2 * L * r) - e(-L * r / 2)) / (1 - e(-2 * L * r))).sum()
    c = -(arb.log(arb.pi()) + arb.const_euler() + arb.log(1 - arb.exp(-4 * L)))
    u = r - 1
    nu_cosh = cheb_moments(u, cosh_w, D)
    Tm1 = [arb(-1) ** j for j in range(D + 1)]
    nu_arch = [a - b for a, b in zip([Tm1[j] * k_w.sum() for j in range(D + 1)], cheb_moments(u, k_w, D))]
    X = int(float((2 * L).exp().upper()))
    pp, lam = [], []
    for p in primes_upto(X):
        q_ = int(p)
        while q_ <= X:
            if arb(q_).log() < 2 * L: pp.append(q_); lam.append(arb(int(p)).log())
            q_ *= int(p)
    if pp:
        xs = np.array([arb(n).log() / L - 1 for n in pp], dtype=object)
        cn = np.array([-2 * l / arb(n).sqrt() for n, l in zip(pp, lam)], dtype=object)
        nu_p = cheb_moments(xs, cn, D)
    else:
        nu_p = [arb(0)] * (D + 1)
    n = D + 1
    def assemble(nu):
        Q = arb_mat(O.K, O.K)
        for i, ui in enumerate(O.u):
            T = cheb_moments(np.array([ui], dtype=object), np.array([arb(1)], dtype=object), D)
            wi = (nu[0] / 2 + sum(nu[j] * T[j] for j in range(1, n))) * 2 / n
            Q += O.G[i] * wi
        return Q
    I = arb_mat(O.K, O.K); [I.__setitem__((i, i), arb(1)) for i in range(O.K)]
    if parts:
        return dict(cosh=assemble(nu_cosh), arch=assemble(nu_arch) + I * d, const=I * c, primes=assemble(nu_p))
    return assemble([a + b + p for a, b, p in zip(nu_cosh, nu_arch, nu_p)]) + I * (d + c)

def lowest(Q, its=6):
    """Lowest two eigenvalues by subspace inverse iteration, with Rayleigh quotients and residual radii."""
    # iterate on midpoints: ball inversion of a matrix with lambda0 ~ 1e-57 inflates radii to nan; the final
    # Rayleigh quotient and residual below use the true Q, so the reported value stays honest
    K = Q.nrows(); mid = lambda M: arb_mat([[x.mid() for x in M.entries()[i * M.ncols():(i + 1) * M.ncols()]] for i in range(M.nrows())])
    I = arb_mat(K, K); [I.__setitem__((i, i), arb(1)) for i in range(K)]
    Qm = mid(Q); Qi = Qm.solve(I, algorithm="approx")  # ball inv() midpoint is inaccurate here (review: off by 50 orders at K=64)
    del Qm, I; Qi = mid(Qi)  # free copies early: at K ~ 3000 and 9000+ bits each K x K copy is ~10-15 GB
    V = arb_mat([[arb(1) / (i + 1) ** j for j in range(2)] for i in range(K)])
    def orth(V):  # Gram-Schmidt on two columns
        a = arb_mat([[V[i, 0]] for i in range(K)]); na = (a.transpose() * a)[0, 0].sqrt(); a = a / na
        b = arb_mat([[V[i, 1]] for i in range(K)]); b = b - a * (a.transpose() * b)[0, 0]; b = b / (b.transpose() * b)[0, 0].sqrt()
        return a, b
    def rq(V):
        out = []
        for v in orth(V):
            v = arb_mat([[x.mid()] for x in v.entries()]); v = v / (v.transpose() * v)[0, 0].sqrt()
            Qv = Q * v; rho = (v.transpose() * Qv)[0, 0]; res = Qv - v * rho
            out.append((rho, (res.transpose() * res)[0, 0].sqrt()))
        return out
    for it in range(1, 400):  # at least `its` steps, then until both residuals are tiny (slow only when lambda1/lambda2 ~ 1)
        a, b = orth(mid(Qi * arb_mat([[V[i, 0], V[i, 1]] for i in range(K)])))
        V = arb_mat([[a[i, 0], b[i, 0]] for i in range(K)])
        if it >= its and it % 5 == 0:
            out = rq(V)
            if all(e.mid() < 1e-6 * abs(r.mid()) for r, e in out): break
    out = rq(V)
    assert all(e.mid() < 1e-6 * abs(r.mid()) for r, e in out), f"inverse iteration failed: {[(r.mid().str(3), e.mid().str(3)) for r, e in out]}"
    return out  # [(lambda0, residual), (lambda1, residual)]

if __name__ == "__main__":
    import sys
    K, par, prec = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); ctx.prec = prec
    t = time.time(); O = Overlap(K, par); print(f"overlap K={K} par={par} prec={prec}: {time.time()-t:.1f}s", flush=True)
    for L in sys.argv[4:]:
        t = time.time(); Q = weil_matrix(O, L); (l0, e0), (l1, e1) = lowest(Q)
        print(f"L={L}: lambda0 = {l0.mid().str(8)} (residual {e0.mid().str(2)}), lambda1 = {l1.mid().str(8)}  {time.time()-t:.1f}s", flush=True)
