# Folding test: phi = lowest odd eigvec, psi = phi on (0,L). Then
#   Q(phi*sgn) - Q(phi) = 2 int_0^{2L} W(x) (psi*psi)(x) dx - sum 4 Lam(n)/sqrt(n) (psi*psi)(log n)  =: Delta*||phi||^2
# lam0_even <= lam0_odd + Delta.  Conjecture F: Delta < 0  (=> even strictly below odd).
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))  # for convergence.py
import mpmath as mp
from convergence import galerkin, legendre_norm
from flint import arb
mp.mp.dps = 60
L, D = float(sys.argv[1]), int(sys.argv[2])
Q = galerkin(L, D)
od = list(range(1, D, 2)); evn = list(range(0, D, 2))
def low(idx):
    A = mp.matrix([[mp.mpf(Q[i, j].mid().str(55, radius=False)) for j in idx] for i in idx])
    e, V = mp.eigsy(A); k = min(range(len(idx)), key=lambda k: e[k])
    return e[k], [V[i, k] for i in range(len(idx))]
lo, v = low(od); le, _ = low(evn)
def psi(x):  # odd eigvec restricted to (0,L), basis e_m(x) = Pbar_m(x/L)/sqrt(L)
    if x <= 0 or x >= L: return mp.mpf(0)
    P = legendre_norm(D, arb(x / L))
    return sum(v[i] * mp.mpf(P[m].mid().str(50, radius=False)) for i, m in enumerate(od)) / mp.sqrt(L)
mp.mp.dps = 30
s = 1 if mp.quad(psi, [0, L / 2]) > 0 else -1
conv = lambda x: mp.quad(lambda u: psi(u) * psi(x - u), [max(0, x - L), min(x, L)]) if 0 < x < 2 * L else mp.mpf(0)
W = lambda x: 4 * mp.cosh(x / 2) - 2 * mp.exp(-x / 2) / (1 - mp.exp(-2 * x))
smooth = 2 * mp.quad(lambda x: W(x) * conv(x), mp.linspace(0, 2 * L, 9))
pp = [(n, p) for n in range(2, int(mp.exp(2 * L)) + 1) for p in [min(d for d in range(2, n + 1) if n % d == 0)] if p ** round(mp.log(n) / mp.log(p)) == n]
primes = sum(4 * mp.log(p) / mp.sqrt(n) * conv(mp.log(n)) for n, p in pp)
minpsi = min(s * psi(L * k / 200) for k in range(1, 200))
print(f"L={L} lam0even={mp.nstr(le,6)} lam0odd={mp.nstr(lo,6)} smooth={mp.nstr(smooth,8)} primes={mp.nstr(-primes,8)} Delta={mp.nstr(smooth-primes,6)} psi_onesigned={minpsi>-1e-12}")
