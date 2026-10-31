"""Per-window parameters for the simple-even grid, interpolated from Elston's three published runs
(1.19: delta 0.3, beta_k 64, T 700, N 680, gl 56; 1.3: 0.3, 80, 1240, 1250, 64; 1.495: 0.2, 120, 2870, 3200, 96).
Omega = T - beta_k/delta is log-interpolated (x1.03); N = 0.8 L (Omega + 64/delta) + 20 as before.
beta_k is the smallest value whose entry radius (~50 eps_side^2, eps_side = 2/(pi beta I0(beta)), fitted to Elston's runs)
times N stays below 1e-6 of the target eigenvalue (log-interpolated odd lambda_0; the binding sector).
First grid run kept beta_k = 64 everywhere: entry radius 2e-56 x N swamped lambda ~ 6e-55 at 1.2785 (certified bound negative)."""
import math
from mpmath import besseli, pi
PTS = [(1.19, 486.667, 64, 56), (1.3, 973.333, 80, 64), (1.495, 2270.0, 120, 96)]
LAM = [(1.25, 2.67e-51), (1.2785, 5.96e-55), (1.3, 5.36e-58), (1.495, 8.26e-92), (1.5, 8.2e-93)]

def interp(L, xs, f=lambda v: v, g=lambda v: v):
    (a, va), (b, vb) = (xs[0], xs[1]) if L <= xs[1][0] else (xs[-2], xs[-1])
    for i in range(len(xs) - 1):
        if xs[i][0] <= L <= xs[i + 1][0]: (a, va), (b, vb) = xs[i], xs[i + 1]
    return g(f(va) + (f(vb) - f(va)) * (L - a) / (b - a))

def params(L):
    p0, p1 = (PTS[0], PTS[1]) if L <= 1.3 else (PTS[1], PTS[2])
    t = (L - p0[0]) / (p1[0] - p0[0])
    Om = math.exp(math.log(p0[1]) + (math.log(p1[1]) - math.log(p0[1])) * t) * 1.03
    gl = max(int(math.ceil(p0[3] + (p1[3] - p0[3]) * t)), 56)
    delta = 0.3 if L <= 1.31 else 0.2
    N = int(0.8 * L * (Om + 64 / delta)) + 20
    lam = interp(L, LAM, math.log, math.exp)
    bk = 64
    while N * 50 * float(2 / (pi * bk * besseli(0, bk))) ** 2 > 1e-6 * lam: bk += 1
    T = int(math.ceil((Om + bk / delta) / 10) * 10)
    return dict(a=f"{L}", delta=str(delta), tsharp=str(T), nmodes=N, gl=gl, prec=400 if L < 1.35 else 600, beta_k=bk)

GRID = [1.2205, 1.25, 1.2785, 1.3055, 1.3315, 1.3565, 1.3805, 1.4035, 1.4255, 1.447, 1.468, 1.4885, 1.5]
if __name__ == "__main__":
    for L in GRID:
        p = params(L)
        print(f"--a {p['a']} --delta {p['delta']} --tsharp {p['tsharp']} --nmodes {p['nmodes']} --gl {p['gl']} --prec {p['prec']} --beta-k {p['beta_k']}")
