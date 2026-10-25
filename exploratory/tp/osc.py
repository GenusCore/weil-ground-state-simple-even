# Oscillation (Sturm/Markov-system) test for the Weil form Q_L: does xi_n have exactly n sign changes?
import os
import sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))  # for convergence.py
import mpmath as mp
from convergence import galerkin, legendre_norm
from flint import arb
mp.mp.dps = 140
L, D, nshow = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
t = time.time()
Q = galerkin(L, D)
A = mp.matrix(D, D)
for i in range(D):
    for j in range(D):
        A[i, j] = mp.mpf(Q[i, j].mid().str(130, radius=False))
ev, V = mp.eigsy(A)
order = sorted(range(D), key=lambda k: ev[k])
xs = [L * (-1 + 2 * (k + 0.5) / 4000) for k in range(4000)]
Ps = [[mp.mpf(p.mid().str(40, radius=False)) for p in legendre_norm(D, arb(x / L))] for x in xs]
for r, k in enumerate(order[:nshow]):
    v = [V[i, k] for i in range(D)]
    par = 'even' if sum(v[i]**2 for i in range(1, D, 2)) < 1e-20 else ('odd' if sum(v[i]**2 for i in range(0, D, 2)) < 1e-20 else 'MIX')
    f = [sum(v[i] * P[i] for i in range(D)) for P in Ps]
    fmax = max(abs(y) for y in f)
    s = [y for y in f if abs(y) > 1e-8 * fmax]  # ignore numerically-zero tails
    sc = sum(1 for a, b in zip(s, s[1:]) if a * b < 0)
    zs = [(round(xs[i],2), int(mp.log10(max(abs(y) for y in f[max(i-20,0):i+20])/fmax))) for i in range(len(f)-1) if f[i]*f[i+1] < 0]
    print(f"L={L} n={r} lam={mp.nstr(ev[k], 6)} {par} signchanges={sc} zeros(x,log10 local amp)={zs}", flush=True)
print('sec', round(time.time() - t))
