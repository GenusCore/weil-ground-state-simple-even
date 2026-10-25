"""Derivative ladder: is the j-th state of the staircase ~ the j-th derivative of the ground state k_L
(Gram-Schmidt'ed against the lower states of the same parity)? Baseline: x^j k_L. Non-rigorous Galerkin."""
import os
import sys
from math import ceil, e, exp, pi as fpi
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))  # for convergence.py
import mpmath as mp
import convergence as cv
mp.mp.dps = 150
J = 6

def deriv(c, L):  # coefficients of f' in the orthonormal basis e_m = sqrt((2m+1)/2L) P_m(x/L)
    D = len(c)
    a = [c[m] * mp.sqrt((2*m+1) / (2*L)) for m in range(D)]  # f = sum a_m P_m(s), s = x/L
    b = [mp.mpf(0)] * D
    for m in range(D):
        b[m] = (2*m+1) * sum(a[j] for j in range(m+1, D, 2)) / L
    return [b[m] / mp.sqrt((2*m+1) / (2*L)) for m in range(D)]

def xmul(c, L):  # coefficients of x f: s P_m = ((m+1) P_{m+1} + m P_{m-1})/(2m+1)
    D = len(c)
    a = [c[m] * mp.sqrt((2*m+1) / (2*L)) for m in range(D)]
    b = [mp.mpf(0)] * (D + 1)
    for m in range(D):
        b[m+1] += a[m] * (m+1) / (2*m+1) * L
        if m: b[m-1] += a[m] * m / (2*m+1) * L
    return [b[m] / mp.sqrt((2*m+1) / (2*L)) for m in range(D)]  # drop top term (tiny)

def dot(u, v): return sum(x*y for x, y in zip(u, v))

def gs(v, basis):
    for b in basis: v = [x - dot(v, b) * y for x, y in zip(v, b)]
    n = mp.sqrt(dot(v, v)); return [x / n for x in v]

def one(Ldec):
    L = mp.mpf(Ldec)
    D = min(2 * ceil(e * float(L) * 2 * fpi * exp(2 * float(L)) / 8) + 60, 170)
    Q = cv.galerkin(float(Ldec), D)
    st = []
    for par in (0, 1):
        idx = list(range(par, D, 2))
        M = mp.matrix([[mp.mpf(Q[i, j].mid().str(140, radius=False)) for j in idx] for i in idx])
        E, V = mp.eigsy(M)
        for k in range(len(idx)):
            v = [mp.mpf(0)] * D
            for i, m in enumerate(idx): v[m] = V[i, k]
            st.append((E[k], par, v))
    st.sort(key=lambda s: s[0])
    k0 = st[0][2]
    rows = []
    for name, op in (("d/dx", deriv), ("x*", xmul)):
        lad, cur, basis = [], k0, {0: [k0], 1: []}
        for j in range(1, J):
            cur = op(cur, L)
            par = j % 2
            w = gs(cur, basis[par]); basis[par].append(w)
            lad.append(w)
        # compare ladder j with staircase state j
        out = []
        for j, w in enumerate(lad, 1):
            lam, par, v = st[j]
            out.append(f"{1 - abs(dot(w, v)):.1e}" if par == j % 2 else "parity!")
        rows.append((name, out))
    return Ldec, D, [mp.nstr(s[0], 3) for s in st[:J]], rows

if __name__ == "__main__":
    with Pool(5) as p:
        for Ldec, D, lams, rows in p.imap(one, ["0.6", "0.8", "1.0", "1.2", "1.4"]):
            print(f"L={Ldec} D={D} staircase {lams}")
            for name, out in rows: print(f"   1-|cos(state j, GS({name})^j k)|  j=1..{J-1}: {out}")
            sys.stdout.flush()
