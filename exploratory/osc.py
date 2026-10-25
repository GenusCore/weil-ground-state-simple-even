"""(1) Oscillation test: does the k-th eigenfunction (both parities merged) have exactly k sign changes in (-L, L)?
(2) Does the ground state converge to Riemann's Phi kernel (Xi(t) = int Phi(u) e^{itu} du)?  Non-rigorous Galerkin."""
import os
import sys
from math import ceil, e, exp, pi as fpi
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))  # for convergence.py
import mpmath as mp
import convergence as cv
mp.mp.dps = 120

def Phi(u):
    return mp.nsum(lambda n: (2 * mp.pi**2 * n**4 * mp.e**(9*u) - 3 * mp.pi * n**2 * mp.e**(5*u)) * mp.e**(-mp.pi * n**2 * mp.e**(4*u)), [1, mp.inf])

def one(Ldec):
    L = float(Ldec)
    D = min(2 * ceil(e * L * 2 * fpi * exp(2 * L) / 8) + 60, 200)
    Q = cv.galerkin(L, D)
    pairs = []
    for par in (0, 1):
        idx = list(range(par, D, 2))
        M = mp.matrix([[mp.mpf(Q[i, j].mid().str(110, radius=False)) for j in idx] for i in idx])
        E, V = mp.eigsy(M)
        for k in range(len(idx)):
            pairs.append((E[k], par, [V[i, k] for i in range(len(idx))], idx))
    pairs.sort(key=lambda p: p[0])
    Lm = mp.mpf(Ldec); xs = [Lm * (-1 + 2 * mp.mpf(i) / 2000) for i in range(1, 2000)]
    def fval(c, idx, x):
        s = x / Lm; P = [mp.mpf(1), s]
        for k in range(1, max(idx)):
            P.append(((2*k+1) * s * P[k] - k * P[k-1]) / (k+1))
        return sum(cm * P[m] * mp.sqrt((2*m+1) / (2*Lm)) for cm, m in zip(c, idx))
    rows = []
    for k, (lam, par, c, idx) in enumerate(pairs[:14]):
        f = [fval(c, idx, x) for x in xs]
        mx = max(abs(v) for v in f)
        sig = [v for v in f if abs(v) > mx * mp.mpf(10)**-25]
        changes = sum(1 for a, b in zip(sig, sig[1:]) if a * b < 0)
        rows.append((k, mp.nstr(lam, 4), "even" if par == 0 else "odd", changes))
        if k == 0:
            g = f
    # ground state vs Phi on the same grid (both L2-normalised on [-L,L]); sign-fix
    ph = [Phi(x) for x in xs]
    h = 2 * Lm / 2000
    n1 = mp.sqrt(sum(v*v for v in g) * h); n2 = mp.sqrt(sum(v*v for v in ph) * h)
    sgn = 1 if g[999] > 0 else -1
    ov = sum(sgn * a * b for a, b in zip(g, ph)) * h / (n1 * n2)
    return Ldec, D, rows, float(ov), float(sgn * g[999] / n1), float(ph[999] / n2)

if __name__ == "__main__":
    with Pool(5) as p:
        for Ldec, D, rows, ov, gpk, ppk in p.imap(one, ["0.3", "0.6", "0.8", "1.0", "1.2"]):
            print(f"L={Ldec} D={D}  overlap(ground, Phi)={ov:.10f}  peak ground={gpk:.5f} Phi={ppk:.5f}")
            print("   k, lambda, parity, sign changes:", rows, flush=True)
