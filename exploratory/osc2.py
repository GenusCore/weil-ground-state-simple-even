"""Where are the sign changes of eigenfunctions k=0..12 at L=1.0, and are they stable under D? (non-rigorous)"""
import os
import sys
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))  # for convergence.py
import mpmath as mp
import convergence as cv
mp.mp.dps = 150

def one(D):
    Ldec = "1.0"; Lm = mp.mpf(Ldec)
    Q = cv.galerkin(1.0, D)
    pairs = []
    for par in (0, 1):
        idx = list(range(par, D, 2))
        M = mp.matrix([[mp.mpf(Q[i, j].mid().str(140, radius=False)) for j in idx] for i in idx])
        E, V = mp.eigsy(M)
        pairs += [(E[k], par, [V[i, k] for i in range(len(idx))], idx) for k in range(len(idx))]
    pairs.sort(key=lambda p: p[0])
    xs = [Lm * mp.mpf(i) / 4000 for i in range(0, 4000)]  # x >= 0 half (parity gives the rest)
    out = []
    for k, (lam, par, c, idx) in enumerate(pairs[:13]):
        def f(x):
            s = x / Lm; P = [mp.mpf(1), s]
            for j in range(1, max(idx)):
                P.append(((2*j+1) * s * P[j] - j * P[j-1]) / (j+1))
            return sum(cm * P[m] * mp.sqrt((2*m+1) / (2*Lm)) for cm, m in zip(c, idx))
        v = [f(x) for x in xs]; mx = max(abs(t) for t in v)
        zs = [(float(xs[i]), float(mp.log10(max(abs(v[i]), abs(v[i+1])) / mx))) for i in range(1, len(v) - 1) if v[i] * v[i+1] < 0]
        out.append((k, mp.nstr(lam, 4), par, [(round(a, 3), round(b, 1)) for a, b in zs]))
    return D, out

if __name__ == "__main__":
    with Pool(2) as p:
        for D, out in p.imap(one, [92, 140]):
            print("D =", D)
            for k, lam, par, zs in out:
                print(f"  k={k:2d} lam={lam:>10} {'even' if par==0 else 'odd '}  zeros x>0 (x, log10|f|/max): {zs}")
            sys.stdout.flush()
