#!/usr/bin/env python3
"""Window-compressed prime-shift norm versus the torus supremum (fail-fast, numerical).

Zhu's envelope bounds the prime comb P(t) = Σ (2c_n/√n) cos(t log n) by its supremum
A = Σ 2|c_n|/√n, which Kronecker's theorem says the line t·(log p)_p approaches on the prime
torus. As a multiplier, P acts in u-space as Σ (c_n/√n)(τ_{log n} + τ_{−log n}). On functions
supported in [−b, b], its top eigenvalue μ(b) = sup ⟨f, P f⟩/‖f‖² is the prime grid's
weighted adjacency spectrum seen through the window, and it can be well below A.

This script estimates μ(b) by piecewise-constant Galerkin (a lower bound that converges as the
cells shrink), using power iteration on P + A·I. It is numerical only; a rigorous upper bound
needs a Schur test.

Usage: .venv/bin/python scripts/grid_norm.py --function zeta --b 0.8,0.9,1.0,1.19 --k 1000,2000,4000
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from progress import Job  # noqa: E402


def mangoldt_weights(n_max):
    out = []
    for n in range(2, n_max + 1):
        m, p = n, None
        for q in range(2, n + 1):
            if m % q == 0:
                p = q
                break
        while m % p == 0:
            m //= p
        if m == 1:
            out.append((math.log(n), math.log(p) / math.sqrt(n)))
    return out


def weights(function, b, comb_a=None):
    """Comb weights for log n < 2b, or for log n < 2·comb_a when the operator is the comb of a
    test function supported in [−comb_a, comb_a] but acts on a wider window [−b, b]."""
    n_max = int(math.floor(math.exp(2 * (comb_a if comb_a is not None else b))))
    if function == "zeta":
        lim = 2 * (comb_a if comb_a is not None else b)
        return [(l, wn) for l, wn in mangoldt_weights(n_max) if l < lim]
    import mpmath as mp
    if function == "dh":
        import connes_letter_mp as cl
        c = cl.dh_log_derivative(n_max)
    elif function == "tstar":
        import conductor5_family_mp as c5
        c = c5.log_derivative(c5.coefficients(c5.tstar(), n_max), n_max)
    elif function == "z1":
        import epstein_connes_mp as ec
        c = ec.log_derivative(ec.coefficients("Z1", n_max), n_max)
    else:
        raise ValueError(function)
    return [(math.log(n), float(c[n]) / math.sqrt(n)) for n in range(2, n_max + 1) if float(c[n]) != 0.0]


def top_eigenvalue(w, b, K, iters=4000, tol=1e-12):
    h = 2 * b / K
    # Toeplitz Galerkin matrix: G[i][j] = Σ_n w_n [tri((j−i)h − ℓ) + tri((j−i)h + ℓ)]/h, tri(d) = max(0, h − |d|).
    diag = {}
    for ell, wn in w:
        for sgn in (1, -1):
            centre = sgn * ell / h
            for d in range(int(math.floor(centre)) - 1, int(math.ceil(centre)) + 2):
                ov = max(0.0, h - abs(d * h - sgn * ell))
                if ov > 0 and abs(d) < K:
                    diag[d] = diag.get(d, 0.0) + wn * ov / h
    A = sum(2 * abs(wn) for _, wn in w)
    offs = sorted(diag.items())
    v = [1.0] * K
    lam = 0.0
    for it in range(iters):
        u = [A * x for x in v]
        for d, g in offs:
            if d >= 0:
                for i in range(K - d):
                    u[i] += g * v[i + d]
            else:
                for i in range(-d, K):
                    u[i] += g * v[i + d]
        nrm = math.sqrt(sum(x * x for x in u))
        new = sum(x * y for x, y in zip(u, v)) / sum(x * x for x in v) - A
        v = [x / nrm for x in u]
        if it > 50 and abs(new - lam) < tol * abs(new):
            lam = new
            break
        lam = new
    return lam, A, it, v


def schur_weight(n_list, b, K, phi, iters=3000):
    """Collatz–Wielandt iteration of the Schur majorant T(φ)_i = Σ_n w_n (max φ over the two cells
    hit by i + ℓ_n/h, plus the same for i − ℓ_n/h). It is monotone and homogeneous, and its fixed
    direction makes the Schur ratio as uniform as piecewise-constant weights allow."""
    h = 2 * b / K
    sh = []
    for n in n_list:
        m, p = n, next(d for d in range(2, n + 1) if n % d == 0)
        while m % p == 0:
            m //= p
        sh.append((int(math.floor(math.log(n) / h)), math.log(p) / math.sqrt(n)))
    v = [max(abs(x), 1e-12) for x in phi]
    for _ in range(iters):
        u = [0.0] * K
        for q, wn in sh:
            for i in range(K):
                a1 = v[i + q] if 0 <= i + q < K else 0.0
                a2 = v[i + q + 1] if 0 <= i + q + 1 < K else 0.0
                b1 = v[i - q - 1] if 0 <= i - q - 1 < K else 0.0
                b2 = v[i - q] if 0 <= i - q < K else 0.0
                u[i] += wn * (max(a1, a2) + max(b1, b2))
        s_ = max(u)
        v = [max(x / s_, 1e-300) for x in u]
    return v


def _exact_max(vals):
    """Largest of a few exact (radius-zero) arb values, compared exactly, without float conversion."""
    from flint import arb
    best = arb(0)
    for v in vals:
        if v > best:
            best = v
    return best


def schur_upper(n_list, b, K, phi):
    """Rigorous upper bound on ‖P‖ on L²[−b, b] (nonnegative weights) by Schur's test with the
    piecewise-constant weight φ on K equal cells: for u in cell i, the shifted cell i ± ℓ/h meets at
    most the two cells i ± q, i ± (q + 1) with q = ⌊ℓ/h⌋, certified in Arb. φ can be any positive
    vector; the Galerkin Perron vector makes the bound tight."""
    from flint import arb, ctx
    old_prec = ctx.prec
    ctx.prec = max(200, old_prec)
    h = 2 * arb(b) / K
    shifts = []
    for n in n_list:
        ell = arb(n).log()
        r = ell / h
        q = int(math.floor(float(r.mid())))
        if not (arb(q) < r < arb(q + 1)):
            raise ValueError(f"cannot certify floor(log {n}/h)")
        m, p = n, next(d for d in range(2, n + 1) if n % d == 0)
        while m % p == 0:
            m //= p
        assert m == 1, "Schur test here is for ζ (prime powers, weights Λ(n)/√n ≥ 0)"
        shifts.append((q, arb(p).log() / arb(n).sqrt()))
    ph = [arb(float(x)) for x in phi]  # floats are exact in arb
    if min(float(x) for x in phi) <= 0:
        raise ValueError("phi must be positive")
    worst = arb(0)  # exact; the row maximum is selected by exact comparison of upper endpoints
    for i in range(K):
        acc = arb(0)
        for q, wn in shifts:
            up = [ph[j] for j in (i + q, i + q + 1) if 0 <= j < K]
            dn = [ph[j] for j in (i - q - 1, i - q) if 0 <= j < K]
            m_up = _exact_max(up)
            m_dn = _exact_max(dn)
            acc += wn * (m_up + m_dn)
        ratio = acc / ph[i]
        u = ratio.upper()  # exact point; comparisons between exact points are certain
        if u > worst:
            worst = u
    ub = worst.upper()
    ctx.prec = old_prec
    return ub


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--function", choices=["zeta", "dh", "tstar", "z1"], default="zeta")
    ap.add_argument("--b", default="0.8,0.9,1.0,1.19")
    ap.add_argument("--k", default="1000,2000,4000")
    ap.add_argument("--schur", action="store_true", help="also certify an upper bound by Schur's test (ζ only)")
    ap.add_argument("--schur-iters", type=int, default=2000)
    ap.add_argument("--comb-a", type=float, default=None, help="restrict the comb to log n < 2·comb_a (window stays [−b, b])")
    args = ap.parse_args()
    bs = [float(s) for s in args.b.split(",")]
    ks = [int(s) for s in args.k.split(",")]
    job = Job(f"grid norm {args.function}", total=len(bs) * len(ks), args=vars(args))
    for b in bs:
        w = weights(args.function, b, args.comb_a)
        T1 = None
        for K in ks:
            mu, A, it, vec = top_eigenvalue(w, b, K)
            extra = ""
            if args.schur and args.function == "zeta":
                n_list = [n for n in range(2, int(math.exp(2 * b)) + 1) if any(abs(math.log(n) - l) < 1e-12 for l, _ in w)]
                ub = schur_upper(n_list, b, K, schur_weight(n_list, b, K, vec, iters=args.schur_iters))
                extra = f"  Schur upper bound (rigorous) = {ub.str(10, radius=False)}"
            job.result(f"b={b} K={K}: mu={mu:.6f}  A={A:.6f}  mu/A={mu / A:.4f}  "
                       f"T1(A)=2πe^A={2 * math.pi * math.exp(A):.1f}  T1(mu)={2 * math.pi * math.exp(mu):.1f}  iters={it}{extra}")
            job.step(f"b={b} K={K}")
    job.done()


if __name__ == "__main__":
    main()
