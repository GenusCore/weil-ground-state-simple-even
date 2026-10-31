"""Weil matrix Q (same as hp.weil_matrix) with no overlap tensor: O(kmax^2) recurrence on the moments.
A_ab = Int_0^2 C_ab(r) dnu(r),  C_ab(r) = Int_{-1}^{1-r} P_a(t) P_b(t+r) dt,  u = r - 1.
Pointwise identity (integrate (2a+1)P_a = (P_{a+1}-P_{a-1})' by parts, a,b >= 1):
    C_{a+1,b} = C_{a-1,b} + (2a+1)/(2b+1) (C_{a,b-1} - C_{a,b+1})
Seeds (exact, in Legendre moments mu_k = Int P_k(u) dnu):  I_k(u) = Int_u^1 P_k,
    C_0b = I_b,  C_a0 = (-1)^a I_a,  C_1b = Int_u^1 s P_b ds - (1+u) I_b.
G_ab = sqrt((2a+1)(2b+1))/2 * (C_ab + C_ba)/2 = nrm_a nrm_b C_ab (a+b even), so Q = nrm A nrm + (d + c) I."""
import sys, time, numpy as np
from flint import arb, arb_mat, ctx
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import hp
from fractions import Fraction

def I(k):  # Int_u^1 P_k as {legendre index: rational coeff as (num, den)}
    return {0: (1, 1), 1: (-1, 1)} if k == 0 else {k - 1: (1, 2 * k + 1), k + 1: (-1, 2 * k + 1)}

def lin(d, mu):
    return sum((mu[k] * p / q for k, (p, q) in d.items()), arb(0))

def mulu(d):  # u P_k = ((k+1) P_{k+1} + k P_{k-1}) / (2k+1), coefficients kept as arb-free fractions
    from fractions import Fraction as F
    out = {}
    for k, (p, q) in d.items():
        c = F(p, q) / (2 * k + 1)
        out[k + 1] = out.get(k + 1, 0) + c * (k + 1)
        if k: out[k - 1] = out.get(k - 1, 0) + c * k
    return {k: (v.numerator, v.denominator) for k, v in out.items() if v}

def add(*ds):
    from fractions import Fraction as F
    out = {}
    for s, d in ds:
        for k, (p, q) in d.items(): out[k] = out.get(k, 0) + s * F(p, q)
    return {k: (v.numerator, v.denominator) for k, v in out.items() if v}

def A_matrix(mu, kmax):
    """A[a][b] for a <= b <= kmax, a + b even (dict of rows).  For a + b even C_ab = C_ba (c_ab(-r) = (-1)^(a+b) c_ab(r)),
    and the recurrence is stable on b >= a (forward errors only grow below the diagonal), so only the upper triangle is built."""
    B = 2 * kmax + 1
    row = {0: np.array([lin(I(b), mu) if b % 2 == 0 else None for b in range(B + 1)], dtype=object)}
    r1 = [None] * (B + 1)
    for b in range(1, B + 1, 2):
        sPb = add((Fraction(b + 1, 2 * b + 1), I(b + 1)), (Fraction(b, 2 * b + 1), I(b - 1)))
        r1[b] = lin(add((1, sPb), (-1, I(b)), (-1, mulu(I(b)))), mu)
    row[1] = np.array(r1, dtype=object)
    inv = np.array([arb(1) / (2 * b + 1) for b in range(B + 2)], dtype=object)
    for a in range(1, kmax):
        n = a + 1; bs = np.arange(n, B - n + 1, 2)  # row n needed on n <= b <= B - n, b = n mod 2
        new = np.empty(B + 1, dtype=object)
        new[bs] = row[a - 1][bs] + (2 * a + 1) * inv[bs] * (row[a][bs - 1] - row[a][bs + 1])
        row[n] = new
        if a - 1 >= 2: row[a - 1] = row[a - 1][: kmax + 1]  # trim, keep memory O(kmax^2)
    return row

def leg_moments(xs, c, N):
    """sum_n c_n P_j(x_n), j = 0..N, with radii that do not blow up.  The ball three-term recurrence inflates radii by
    (1+sqrt 2)^j (majorant of the error recurrence), ~1.27 bits per degree; instead run the Chebyshev recurrence on exact
    midpoints, re-centring every step, and add a rigorous bound: the error of T_j = 2x T_{j-1} - T_{j-2} with per-step
    rounding delta is sum_i delta_i U_{j-i-1}(x), |U_m| <= m+1, so <= j^2 delta / 2; node radius adds j^2 rad(x) (Markov).
    Then P_n = sum_m lam_m lam_{n-m} (2 - [n=2m]) T_{n-2m}, lam_m = binom(2m,m)/4^m: nonnegative, row sums 1, no amplification."""
    xm = np.array([arb(x.mid()) for x in xs], dtype=object)
    xr = max(x.rad() for x in xs); cabs = sum(abs(v).upper() for v in c)
    delta = arb(2) ** (8 - ctx.prec)  # bound on per-step rounding of 2 x T1 - T0 (|x|, |T| <= 1)
    unit = arb(0, 1); err = lambda j: unit * (cabs * j * j * (delta / 2 + xr)).upper()
    T0 = np.ones(len(xs), dtype=object) * arb(1); T1 = xm.copy()
    t = [c.sum(), (c * T1).sum() + err(1)]
    for j in range(2, N + 1):
        nxt = 2 * xm * T1 - T0
        assert max(v.rad() for v in nxt) <= delta, "per-step rounding exceeds the assumed bound delta"
        T0, T1 = T1, np.array([arb(v.mid()) for v in nxt], dtype=object)
        t.append((c * T1).sum() + err(j))
    lam = [arb(1)]
    for m in range(1, N + 1): lam.append(lam[-1] * (2 * m - 1) / (2 * m))
    return [sum((lam[m] * lam[n - m] * (1 if 2 * m == n else 2) * t[n - 2 * m] for m in range(n // 2 + 1)), arb(0)) for n in range(N + 1)]

def quad_bounds(L, n, N):
    """Rigorous Gauss-Legendre error bounds for the smooth moments (u in [-1, 1], r = u + 1).
    Trefethen ATAP Thm 19.3: if f is analytic in the Bernstein ellipse E_rho with |f| <= M, the n-point Gauss error is
    <= (64/15) M rho^(-2n+2) / (rho^2 - 1) (n-point rule; ATAP states rho^(-2n) for n+1 points, see gl_correction.py).  On E_rho, |P_j| <= rho^j (Laplace integral).  The weights are
    cosh: 4L cosh(Lr/2) (entire); arch: h(r) * (P_j(-1) - P_j(u)) / r with h = r k(r) = 2L r e^(-Lr/2) / (1 - e^(-2Lr)),
    whose poles are at r = i pi m / L, m != 0; diagonal d: 2L (e^(-2Lr) - e^(-Lr/2)) / (1 - e^(-2Lr)), same poles.
    The maxima on E_rho are enclosed by evaluating on acb balls that cover the ellipse.  The analyticity check
    (poles outside E_rho) is exact.  Returns (eps_j for j = 0..N, eps_d)."""
    from flint import acb
    L = arb(L); pi = arb.pi()
    def inside(rho, m):  # is u = -1 + i pi m / L strictly inside E_rho?
        a = (rho + 1 / rho) / 2; b = (rho - 1 / rho) / 2
        return (1 / a) ** 2 + (pi * m / (L * b)) ** 2 < 1
    rho = arb(1)
    for t in range(1, 400):  # largest rho on a grid that keeps the first pole well outside (score >= 1.2)
        cand = 1 + arb(t) / 200; a = (cand + 1 / cand) / 2; b = (cand - 1 / cand) / 2
        if (1 / a) ** 2 + (pi / (L * b)) ** 2 >= arb("1.2"): rho = cand
        else: break
    assert rho > 1 and not inside(rho, 1)  # higher m poles lie further out
    a = (rho + 1 / rho) / 2; b = (rho - 1 / rho) / 2
    S = 2000; Mc = Mh = Md = arb(0); dmin = None
    for s in range(S):  # cover the ellipse u = a cos t + i b sin t by S balls (|u'(t)| <= a)
        t = 2 * pi * (arb(s) + arb("0.5")) / S
        u = acb(a * t.cos(), b * t.sin()); u = acb(u.real.mid(), u.imag.mid()) + acb(0, 0) ; rad = a * pi / S * arb("1.01")
        uu = acb(arb(u.real.mid(), rad.upper()), arb(u.imag.mid(), rad.upper())); r = uu + 1
        Mc = max(Mc, abs(4 * L * (L * r / 2).cosh()).upper(), key=lambda x: float(x))
        e1, e2 = (-L * r / 2).exp(), (-2 * L * r).exp()
        Mh = max(Mh, abs(2 * L * r * e1 / (1 - e2)).upper(), key=lambda x: float(x))
        Md = max(Md, abs(2 * L * (e2 - e1) / (1 - e2)).upper(), key=lambda x: float(x))
        dd = abs(r).lower(); dmin = dd if dmin is None or float(dd) < float(dmin) else dmin
    assert all(x.is_finite() for x in (Mc, Mh, Md)) and dmin > 0
    C = arb(64) / 15 * rho ** (-2 * n + 2) / (rho ** 2 - 1)  # n-point rule: rho^(-2n+2) (gl_correction.py); was rho^(-2n)
    eps = [(C * (Mc * rho ** j + Mh * (rho ** j + 1) / dmin)).upper() for j in range(N + 1)]
    return eps, (C * Md).upper()

def measure(L, D, Mq=None, rigorous=True, primes=True, pole=True):
    """Legendre moments (in u = r - 1) of nu_L and the diagonal constant d + c; identical quadrature to hp.weil_matrix."""
    L = arb(L); Mq = Mq or 2 * D + 600; N = D + 2
    q = hp.gl(Mq); r = np.array([x + 1 for x, _ in q], dtype=object); w = np.array([y for _, y in q], dtype=object)
    e = lambda z: np.array([arb.exp(v) for v in z], dtype=object)
    cosh_w = w * 4 * L * (e(L * r / 2) + e(-L * r / 2)) / 2
    k_w = w * 2 * L * e(-L * r / 2) / (1 - e(-2 * L * r))
    d = (w * 2 * L * (e(-2 * L * r) - e(-L * r / 2)) / (1 - e(-2 * L * r))).sum()
    c = -(arb.log(arb.pi()) + arb.const_euler() + arb.log(1 - arb.exp(-4 * L)))
    u = r - 1
    mu = leg_moments(u, cosh_w, N) if pole else [arb(0)] * (N + 1)  # pole=False: entire L-functions (dh/)
    ks = k_w.sum(); mk = leg_moments(u, k_w, N)
    mu = [m + (-1) ** j * ks - x for j, (m, x) in enumerate(zip(mu, mk))]  # arch: Int g(0) - g(r) against k_w
    if rigorous:  # add the Gauss-Legendre truncation error (not otherwise in the radii)
        eps, eps_d = quad_bounds(L, Mq, N)
        mu = [m + arb(0, e) for m, e in zip(mu, eps)]; d = d + arb(0, eps_d)
    X = int(float((2 * L).exp().upper())); pp, lam = [], []
    for p in hp.primes_upto(X):
        q_ = int(p)
        while q_ <= X:
            if arb(q_).log() < 2 * L and (primes is True or (primes and int(p) in primes)): pp.append(q_); lam.append(arb(int(p)).log())
            q_ *= int(p)
    if pp and primes:  # primes=False: archimedean form only (explore/shape.py); a set of primes: only those (explore/prime_order.py)
        xs = np.array([arb(n).log() / L - 1 for n in pp], dtype=object)
        cn = np.array([-2 * l / arb(n).sqrt() for n, l in zip(pp, lam)], dtype=object)
        mu = [a + b for a, b in zip(mu, leg_moments(xs, cn, N))]
    return mu, d + c

def weil_matrix(K, par, L, Mq=None, meas=None, primes=True):
    kmax = 2 * K - 2 + par; D = 2 * kmax + 1
    mu, dc = meas or measure(L, D, Mq, primes=primes)
    row = A_matrix(mu, kmax)
    ks = list(range(par, kmax + 1, 2)); nrm = [((2 * k + 1) * arb(1) / 2).sqrt() for k in ks]
    Q = arb_mat(K, K)
    for i, a in enumerate(ks):
        for j, b in enumerate(ks):
            Q[i, j] = nrm[i] * nrm[j] * row[min(a, b)][max(a, b)] + (dc if i == j else 0)
    return Q

if __name__ == "__main__":
    K, par, prec = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); ctx.prec = prec
    for L in sys.argv[4:]:
        t = time.time(); Q = weil_matrix(K, par, L); tb = time.time() - t
        rad = max(float(x.rad()) for x in Q.entries())
        (l0, e0), (l1, e1) = hp.lowest(Q)
        print(f"K={K} par={par} prec={prec} L={L}: build {tb:.2f}s maxrad {rad:.2e} lambda0 = {l0.mid().str(12)} (res {e0.mid().str(2)}) lambda1 = {l1.mid().str(8)}", flush=True)
