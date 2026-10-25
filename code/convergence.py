"""Steps 2-3 of REDUCTION.md (internal note, not released): how close is Connes's explicit k_L to the Weil ground state xi_L?

    rho(L)^2 = (Q(k) - lam0 ||k||^2) / ((lam1 - lam0) ||k||^2)   (Lemma 1: ||k - proj_xi k|| <= rho ||k||)

k_L(x) = e^{x/2} sum_n h(n e^x) on [-L, L], with h the vanishing-integral combination of the prolate functions
h_0, h_4 on [-lam, lam], lam = e^L, bandwidth c = 2 pi lam^2 (Connes, arXiv:2602.04022, eq. (17)).

  step 3 (galerkin): high-precision position-side Galerkin matrix of Q in the Legendre basis (non-rigorous).
  step 2 (certified, L <= 0.8): exact Rayleigh quotient of the even Legendre projection of k_L (certify.rayleigh_upper)
          plus the certified lam0, lam1 bounds in theorem_bounds.json.
"""
import json
import sys
from math import ceil, e as E_
import numpy as np
from flint import arb, arb_mat, ctx

ctx.prec = 512
pi = arb.pi()


def legendre_norm(K, s):  # [Pbar_0(s), ..., Pbar_{K-1}(s)], orthonormal on [-1, 1]
    P = [arb(1), s]
    for k in range(1, K - 1):
        P.append(((2 * k + 1) * s * P[k] - k * P[k - 1]) / (k + 1))
    return [P[k] * ((2 * k + 1) / arb(2)).sqrt() for k in range(K)]


def prolate_even(c, K):
    """Even prolate functions h_0, h_2, h_4 on [-1,1] (bandwidth c) as coefficient vectors in Pbar_0, Pbar_2, ..."""
    ks = list(range(0, 2 * K, 2))
    c2 = arb(c) ** 2
    T = [[arb(0)] * K for _ in range(K)]
    for i, k in enumerate(ks):
        T[i][i] = k * (k + 1) + c2 * arb(2 * k * k + 2 * k - 1) / ((2 * k - 1) * (2 * k + 3))
        if i + 1 < K:
            T[i][i + 1] = T[i + 1][i] = c2 * (k + 1) * (k + 2) / ((2 * k + 3) * ((2 * k + 1) * arb(2 * k + 5)).sqrt())
    Tf = np.array([[float(x) for x in r] for r in T])
    ev = np.sort(np.linalg.eigvalsh(Tf))
    out = []
    for target in ev[:3]:  # inverse iteration with a slightly-off shift
        sh = arb(target) + arb(target) * arb("1e-12") + arb("1e-12")
        M = arb_mat([[T[i][j] - (sh if i == j else 0) for j in range(K)] for i in range(K)])
        x = arb_mat([[1] for _ in range(K)])
        for _ in range(8):
            x = M.solve(x)
            n = sum(x[i, 0] ** 2 for i in range(K)).sqrt()
            x = arb_mat([[(x[i, 0] / n).mid()] for i in range(K)])
        out.append([x[i, 0] for i in range(K)])
    return out  # h_0, h_2, h_4


class KL:
    """Connes's k_L on [-L, L]."""

    def __init__(self, L):
        self.L = arb(L)
        self.lam = self.L.exp()
        self.c = 2 * pi * self.lam ** 2
        self.K = int(float(self.c.mid())) // 2 + 60
        h0, _, h4 = prolate_even(self.c, self.K)
        # vanishing integral: int Pbar_0 = sqrt2, all higher Pbar_k integrate to 0
        self.coef = [h4[0] * a - h0[0] * b for a, b in zip(h0, h4)]
        self.nmax = int(float(self.lam.mid()) ** 2) + 1

    def h(self, y):  # y in [-lam, lam]
        s = y / self.lam
        Pb = legendre_norm(2 * self.K, s)
        return sum(cf * Pb[2 * i] for i, cf in enumerate(self.coef))

    def __call__(self, x):
        u = x.exp()
        return u.sqrt() * sum((self.h(n * u) for n in range(1, self.nmax + 1) if n * u <= self.lam), arb(0))

    def breakpoints(self):  # x where n e^x = lam, inside (-L, L)
        L = self.L
        return sorted([L - arb(n).log() for n in range(2, self.nmax + 1) if L - arb(n).log() > -L])


def project(L, D, k=None, nodes=None):
    """Coefficients of k_L in the orthonormal Legendre basis e_m(x) = sqrt((2m+1)/2L) Pbar-normalised P_m(x/L)."""
    k = k or KL(L)
    Lb = arb(L)
    pts = [-Lb] + k.breakpoints() + [Lb]
    nodes = nodes or D + 40
    gauss = [arb.legendre_p_root(nodes, i, weight=True) for i in range(nodes)]
    a = [arb(0)] * D
    norm2 = arb(0)
    for lo, hi in zip(pts, pts[1:]):
        hw, mid = (hi - lo) / 2, (hi + lo) / 2
        for y, w in gauss:
            x = mid + hw * y
            kx = k(x)
            P = legendre_norm(D, x / Lb)
            for m in range(D):
                a[m] += hw * w * kx * P[m] / Lb.sqrt()
            norm2 += hw * w * kx * kx
    return a, norm2


def galerkin(L, D, nout=None):
    """Position-side Galerkin matrix of Q on span{e_0..e_{D-1}} on [-L, L] (high precision, non-rigorous)."""
    Lb = arb(L)
    M = D + 1  # inner Gauss: exact for polynomial integrands of degree <= 2D
    nout = nout or D + 80
    inner = [arb.legendre_p_root(M, i, weight=True) for i in range(M)]
    outer = [arb.legendre_p_root(nout, i, weight=True) for i in range(nout)]

    def G(x):  # G[j,k] = int e_j(y) e_k(y + x) dy, then symmetrised
        hw = (2 * Lb - x) / 2
        ys = [-Lb + hw * (t + 1) for t, _ in inner]
        E1 = arb_mat([[hw * w * p / Lb for (t, w), p in zip(inner, col)] for col in
                      zip(*[legendre_norm(D, y / Lb) for y in ys])])
        E2 = arb_mat(list(map(list, zip(*[legendre_norm(D, (y + x) / Lb) for y in ys]))))
        Gm = E1 * E2.transpose()
        return (Gm + Gm.transpose()) * arb(0.5)

    G0 = G(arb(0))
    gamma = arb.const_euler()
    Q = G0 * (-(pi.log() + gamma + (1 - (-4 * Lb).exp()).log()))
    w0 = arb(0)
    for t, w in outer:
        x = Lb * (t + 1)  # x in (0, 2L), dx = L dt
        q = (-2 * x).exp()
        Q += G(x) * (Lb * w * (4 * (x / 2).cosh() - 2 * (-x / 2).exp() / (1 - q)))
        w0 += Lb * w * 2 * q / (1 - q)
    Q += G0 * w0
    n = 2
    while arb(n).log() < 2 * Lb:
        p = next(d for d in range(2, n + 1) if n % d == 0)
        m = n
        while m % p == 0:
            m //= p
        if m == 1:
            Q -= G(arb(n).log()) * (2 * arb(p).log() / arb(n).sqrt())
        n += 1
    return Q


def eig_low(Q, idx, k=2):
    """Lowest k Galerkin eigenvalues of Q restricted to the basis indices idx (deflated inverse iteration)."""
    n = len(idx)
    A = [[Q[i, j].mid() for j in idx] for i in idx]
    vals, vecs = [], []
    for _ in range(k):
        B = [[A[i][j] + sum(v[i] * v[j] for v in vecs) for j in range(n)] for i in range(n)]
        Bm = arb_mat(B)
        x = arb_mat([[1] for _ in range(n)])
        for _ in range(80):
            x = Bm.solve(x)
            nr = sum(x[i, 0] ** 2 for i in range(n)).sqrt()
            x = arb_mat([[(x[i, 0] / nr).mid()] for i in range(n)])
        v = [x[i, 0] for i in range(n)]
        vals.append(sum(v[i] * A[i][j] * v[j] for i in range(n) for j in range(n)))
        vecs.append(v)
    return vals, vecs


def step3(L, D):
    Q = galerkin(L, D)
    ev, (v0, _) = eig_low(Q, list(range(0, D, 2)), 2)
    od, _ = eig_low(Q, list(range(1, D, 2)), 1)
    a, nk2 = project(L, D)
    Qk = sum(a[i] * Q[i, j] * a[j] for i in range(D) for j in range(D))
    na2 = sum(x * x for x in a)
    lam0, lam1 = ev[0], min(ev[1], od[0])
    rho2 = (Qk / na2 - lam0) / (lam1 - lam0)
    odd_frac = (sum(a[i] ** 2 for i in range(1, D, 2)) / na2).sqrt()
    cosang = abs(sum(a[2 * i] * v0[i] for i in range(len(v0)))) / na2.sqrt()
    return dict(lam0=lam0, lam1=lam1, Qk=Qk / na2, rho=rho2.sqrt() if rho2 > 0 else rho2, odd=odd_frac,
                sin_angle=(1 - cosang ** 2).sqrt(), proj_loss=(1 - na2 / nk2).sqrt() if nk2 > na2 else arb(0))


def step2(Ldec):
    """Certified rho bound for the even Legendre projection of k_L (L <= 0.8, bounds from theorem_bounds.json)."""
    import certify
    lo, U, l1, od = [arb(s) for s in json.load(open("theorem_bounds.json"))[Ldec]]
    W = certify.Window(Ldec, 0)
    with ctx.workprec(768):
        a, nk2 = project(float(Ldec), 2 * W.N)
        a_even = [a[2 * m].mid() for m in range(W.N)]
        Uk = certify.rayleigh_upper(W, a_even)  # rigorous upper bound on Q(k~)/||k~||^2
        rho2 = (Uk - lo) / (min(l1, od) - U)
        na2 = sum(x * x for x in a_even)
        odd = (sum(a[2 * m + 1] ** 2 for m in range(W.N)) / na2).sqrt()
        Lb = arb(Ldec)
        fourier = (2 * Lb).sqrt() * (Lb / 2).exp() * rho2.sqrt()  # Lemma 2 on |Im z| <= 1/2, per unit ||k~||
    return Uk, rho2.sqrt(), fourier, odd, lo, U


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "step3":
        print(" L     D    lam0          lam1          Q(k)/|k|^2    rho           sin(k,xi0)    odd part")
        for arg in sys.argv[2:]:  # "L" or "L:D1,D2,..." (explicit basis sizes for a convergence check)
            Ld, _, Ds = arg.partition(":")
            L = float(Ld)
            D = 2 * ceil(E_ * L * 2 * np.pi * np.exp(2 * L) / 4 / 2) + 40  # ~ e*L*T*/2 modes + margin
            for DD in ([int(d) for d in Ds.split(",")] if Ds else (D, D + 20)):
                r = step3(L, DD)
                print(f"{Ld:4} {DD:4}  {r['lam0'].str(4):>12}  {r['lam1'].str(4):>12}  {r['Qk'].str(4):>12}"
                      f"  {r['rho'].str(3):>12}  {r['sin_angle'].str(3):>12}  {r['odd'].str(3):>10}", flush=True)
    else:
        print(" L     Q(k~)/|k~|^2 <=   lam0 in            rho <=        |a xi^ - k~^|/|k~| <= (|Im z|<=1/2)   odd part of k_L")
        for Ld in sys.argv[2:]:
            Uk, rho, fo, odd, lo, U = step2(Ld)
            print(f"{Ld:4}  {Uk.str(4):>14}  [{lo.str(3)}, {U.str(3)}]  {rho.upper().str(3):>12}  {fo.upper().str(3):>14}"
                  f"   {odd.str(3)}", flush=True)
