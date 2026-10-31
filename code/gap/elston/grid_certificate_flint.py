#!/usr/bin/env python3
"""Zhu's finite reduction with the window-compressed prime norm (prototype, Arb/FLINT).

Zhu's reduction (arXiv:2608.24827, Thm 1.1; audited in docs/AUDIT_ZHU.md) needs
Ψ(t) = H(t) − P(t) ≥ β on [T♯, ∞), so T♯ ≳ 2π e^{A} with A = sup P (Kronecker returns on the
prime torus). Refinement: split each f supported in [−a, a] with m(t) = 1 − K(t), where K is the
Fourier transform of a bump on [−δ, δ] normalised to K(0) = 1:

  h = f − f * bump,  ĥ = mF,  supp h ⊂ [−b, b], b = a + δ,
  ⟨h, P h⟩ ≤ μ(b) ‖h‖²,       μ(b) = top eigenvalue of the prime shifts compressed to [−b, b],
  Q(f) ≥ pole + (1/π) ∫₀^∞ Φ |F|²,   Φ = Ψ·w + (H − μ)·m²,   w = 1 − m².

For large t, Φ = H − μ − w(P − μ) ≥ H − μ − |w|(A + μ), so the threshold becomes T♯ ≳ 2π e^{μ}.
μ(b) ≈ A/3 numerically (scripts/grid_norm.py).

This prototype assembles the reduced block R' (Zhu's matrix with Ψ → Φ) at midpoint values and
reports its smallest eigenvalue. --mode zhu reproduces the audited R block (Φ = Ψ).
NOT A CERTIFICATE: μ is a Galerkin estimate, the tail envelope is scanned, not proved, and the
eigenvalue is computed in floating point.

Usage:
  .venv/bin/python scripts/grid_certificate_flint.py --mode zhu --a 0.8 --tsharp 200 --nmodes 200
  .venv/bin/python scripts/grid_certificate_flint.py --mode grid --a 0.8 --delta 0.1 --mu 1.62 --tsharp 60 --nmodes 60
"""

import argparse
import math
import os
import sys
import time

import mpmath as mp
from flint import acb, arb, arb_mat, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from progress import Job  # noqa: E402


def comb_terms(a):
    out = []
    for n in range(2, int(math.exp(2 * a)) + 1):
        if math.log(n) >= 2 * a:
            continue
        m, p = n, next(q for q in range(2, n + 1) if n % q == 0)
        while m % p == 0:
            m //= p
        if m == 1:
            out.append((arb(n).log(), 2 * arb(p).log() / arb(n).sqrt()))
    return out


def make_funcs(a, comb, mode, delta, k, mu):
    log_pi = arb.pi().log()
    kk = arb(k) + arb(3) / 2

    def H(t):
        return acb(arb(1) / 4, t / 2).digamma().real - log_pi

    def P(t):
        return sum((c * (t * l).cos() for l, c in comb), arb(0))

    def K(t):
        if mode == "zhu":
            return arb(1)
        x = delta * t
        if x == 0:
            return arb(1)
        nu = arb(k) + arb(1) / 2
        return kk.gamma() * (2 / x) ** nu * acb(x).bessel_j(acb(nu)).real

    def Phi(t):
        h, p = H(t), P(t)
        if mode == "zhu":
            return h - p
        m = 1 - K(t)
        w = 1 - m * m
        return (h - p) * w + (h - mu) * m * m

    return H, P, K, Phi


def sph_j_column(L, t, degrees, extra_bits=200):
    """S_n(t) = √(2L(2n+1)) j_n(L t) by Miller's downward recurrence in plain MPFR (gmpy2),
    normalised with Σ_n (2n+1) j_n(x)² = 1. Numerical (no enclosure): ball arithmetic overestimates
    the recurrence's error growth, so the midpoint is computed at prec + extra_bits."""
    import gmpy2
    nmax = degrees[-1]
    xa = L * t
    if xa == 0:
        vals = {0: 1}
        out = []
        for n in degrees:
            out.append((2 * L * (2 * n + 1)).sqrt() * (arb(1) if n == 0 else arb(0)))
        return out
    with gmpy2.context(gmpy2.get_context(), precision=ctx.prec + extra_bits):
        x = gmpy2.mpfr(xa.mid().str(int((ctx.prec + extra_bits) * 0.302) + 5, radius=False))
        xf = float(x)
        n0 = int(max(nmax, xf)) + 60 + int(4 * xf ** (1 / 3))
        jn1, jn = gmpy2.mpfr(0), gmpy2.mpfr(1)
        raw = [None] * (n0 + 1)
        raw[n0] = jn
        norm = gmpy2.mpfr(2 * n0 + 1)
        for n in range(n0, 0, -1):
            jm = (2 * n + 1) / x * jn - jn1
            jn1, jn = jn, jm
            raw[n - 1] = jm
            norm += (2 * n - 1) * jm * jm
        scale = 1 / gmpy2.sqrt(norm)
        if (gmpy2.sin(x) / x) * raw[0] < 0:
            scale = -scale
        out = []
        for n in degrees:
            mant, ex = gmpy2.mpfr(scale * raw[n]).as_mantissa_exp()
            out.append((2 * L * (2 * n + 1)).sqrt() * (arb(int(mant)) * arb(2) ** int(ex)))
        return out


def smallest_two(Am, job, iters=40, count=2):
    """The `count` (1 or 2) smallest eigenvalues of a symmetric positive (near-singular) matrix by inverse
    iteration with deflation, using FLINT's LU solve (numerical: midpoints only). Returns a list padded with
    None to length 2."""
    N = Am.nrows()
    out, vecs = [], []
    for _ in range(count):
        v = arb_mat([[arb(1) / (1 + i)] for i in range(N)])
        lam = None
        for it in range(iters):
            for u in vecs:
                v = v - u * (u.transpose() * v)[0, 0]
            v = Am.solve(v, nonstop=True)
            v = arb_mat([[v[i, 0].mid()] for i in range(N)])
            for u in vecs:
                v = v - u * (u.transpose() * v)[0, 0]
            nrm = sum((v[i, 0] ** 2 for i in range(N)), arb(0)).sqrt()
            v = arb_mat([[(v[i, 0] / nrm).mid()] for i in range(N)])
            new = (v.transpose() * (Am * v))[0, 0].mid()
            if lam is not None and abs(float(((new - lam) / new).mid())) < 1e-25:
                lam = new
                break
            lam = new
        out.append(mp.mpf(lam.str(40, radius=False)))
        vecs.append(v)
        job.log(f"eigenvalue {len(out)}: {lam.str(15, radius=False)} after {it + 1} iterations")
    return out + [None] * (2 - len(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", choices=["zhu", "grid", "kaiser"], default="kaiser")
    ap.add_argument("--omega", default=None, help="kaiser: low-pass cutoff Ω (transition [Ω − β_K/δ, Ω + β_K/δ])")
    ap.add_argument("--beta-k", default="40", help="kaiser: window parameter β_K (ripple ≈ e^{-β_K})")
    ap.add_argument("--sector", choices=["even", "odd"], default="even")
    ap.add_argument("--a", default="0.8")
    ap.add_argument("--delta", default="0.1")
    ap.add_argument("--k", type=int, default=8, help="bump (1 − u²/δ²)^k")
    ap.add_argument("--mu", default=None, help="compressed norm μ(a + δ) (grid mode)")
    ap.add_argument("--tsharp", type=int, default=200)
    ap.add_argument("--panel", default="0.25")
    ap.add_argument("--gl", type=int, default=32)
    ap.add_argument("--nmodes", type=int, default=200)
    ap.add_argument("--prec", type=int, default=400)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chunk", type=int, default=2048)
    ap.add_argument("--scan", type=int, default=20, help="scan Φ on [T♯, scan·T♯] to set β")
    ap.add_argument("--blocks", default="", help="also report λ₁ of the leading k×k blocks, e.g. 400,500")
    args = ap.parse_args()
    ctx.prec = args.prec
    ctx.threads = args.threads
    mp.mp.prec = args.prec
    a, delta = arb(args.a), arb(args.delta)
    L = a
    comb = comb_terms(float(args.a))
    A_sup = sum((c for _, c in comb), arb(0))
    mu = arb(args.mu) if args.mu else A_sup
    H, P, K, Phi = make_funcs(a, comb, "grid" if args.mode == "kaiser" else args.mode, delta, args.k, mu)
    job = Job(f"grid certificate {args.mode} {args.sector} a={args.a} T{args.tsharp} N{args.nmodes}", args=vars(args))
    job.log(f"prime powers: {[int(round(float(l.exp().mid()))) for l, _ in comb]}; A = {A_sup.str(10)}; mu = {mu.str(10)}")

    # β: Zhu's β* (mode zhu); in grid mode the analytic envelope log(T/2π) − 1/T − μ − (A+μ)·sup|w|,
    # with sup|w| on [T♯, scan·T♯] taken from a scan (prototype).
    T = arb(args.tsharp)
    if args.mode == "zhu":
        beta = (T / (2 * arb.pi())).log() - 1 / T - A_sup
        wsup = arb(0)
    elif args.mode == "kaiser":
        bk = arb(args.beta_k)
        half = bk / delta
        omega = arb(args.omega)
        assert abs(float((omega + half - T).mid())) < 1e-9, "set --tsharp = Ω + β_K/δ"
        wsup = arb(0)  # prototype: ripple ≈ 1/I0(β_K) beyond T♯ neglected
        beta = (T / (2 * arb.pi())).log() - 1 / T - mu
    else:
        wsup = 0.0
        steps = 4000
        for i in range(steps + 1):
            t = T * (1 + (args.scan - 1) * arb(i) / steps)
            m = 1 - K(t)
            wsup = max(wsup, float(abs(1 - m * m).upper()))
        wsup = arb(wsup)
        beta = (T / (2 * arb.pi())).log() - 1 / T - mu - (A_sup + mu) * wsup
    phi_min = None
    for i in range(2001 if args.mode != "kaiser" else 0):
        t = T * (1 + (args.scan - 1) * arb(i) / 2000)
        v = float(Phi(t).lower())
        phi_min = v if phi_min is None else min(phi_min, v)
    job.log(f"beta = {beta.str(12)}; sup|w| on scan = {wsup.str(4)}; min Φ on scan = {phi_min}")
    if beta <= 0:
        job.warn("beta ≤ 0: T♯ too small for the envelope")

    first = 0 if args.sector == "even" else 1
    degrees = [first + 2 * i for i in range(args.nmodes)]
    sign = 1 if args.sector == "even" else -1
    N = len(degrees)

    # Gauss–Legendre nodes on panels of [0, T♯]
    gl = [arb.legendre_p_root(args.gl, i, weight=True) for i in range(args.gl)]
    h = arb(args.panel)
    panels = int(round(args.tsharp / float(args.panel)))
    nodes = []
    for p in range(panels):
        c = h * p + h / 2
        for x, wt in gl:
            nodes.append((c + x * h / 2, wt * h / 2))
    job.total = (len(nodes) + args.chunk - 1) // args.chunk
    if args.mode == "kaiser":
        i0 = acb(bk).bessel_i(0).real

        def W_hat(sv):  # Kaiser spectrum / (2π), entire in s
            z = bk * bk - delta * delta * sv * sv
            r = acb(z).sqrt()
            val = (r.sinh() / r).real if z != 0 else arb(1)
            return 2 * delta / i0 * val / (2 * arb.pi())

        lo, hi = omega - half, omega + half
        band = sorted({float(t.mid()) for t, _ in nodes if lo < t < hi})
        mtab, G, prev = {}, arb(0), float((lo - omega).mid())
        g4 = [arb.legendre_p_root(4, i, weight=True) for i in range(4)]
        for tv in band:
            tau = tv - float(omega.mid())
            aa, bb = arb(prev), arb(tau)
            G += sum((wq * W_hat((aa + bb) / 2 + xq * (bb - aa) / 2) for xq, wq in g4), arb(0)) * (bb - aa) / 2
            mtab[tv] = G
            prev = tau
        job.log(f"kaiser: Ω = {omega.str(6)}, band [{lo.str(6)}, {hi.str(6)}], {len(band)} band nodes, G(end) = {G.str(8)}, 1/I0(β_K) = {(1 / i0).str(3)}")

        def Phi(t, _H=H, _P=P):
            tv = float(t.mid())
            m = arb(0) if tv <= float(lo.mid()) else (arb(1) if tv >= float(hi.mid()) else mtab[tv])
            hh, pp = _H(t), _P(t)
            w = 1 - m * m
            return (hh - pp) * w + (hh - mu) * m * m
    t0 = time.time()
    G = arb_mat(N, N)
    for s in range(0, len(nodes), args.chunk):
        chunk = nodes[s:s + args.chunk]
        cols, scaled = [], []
        for t, wt in chunk:
            col = sph_j_column(L, t, degrees)
            cq = wt * (Phi(t) - beta) / arb.pi()
            cols.append(col)
            scaled.append([cq * v for v in col])
        S = arb_mat([[cols[j][i] for j in range(len(chunk))] for i in range(N)])
        D = arb_mat([[scaled[j][i] for j in range(len(chunk))] for i in range(N)])
        G += D * S.transpose()
        job.step(f"{s + len(chunk)}/{len(nodes)} nodes")
    job.log(f"Gram assembled in {time.time() - t0:.0f}s")

    # pole vector, sign-adjusted
    pole = []
    for n in degrees:
        sg = (-1) ** (n // 2) if args.sector == "even" else (-1) ** ((n - 1) // 2)
        x = L / 2
        i_n = (arb.pi() / (2 * x)).sqrt() * acb(x).bessel_i(acb(n + arb(1) / 2)).real
        pole.append(sg * (2 * L * (2 * n + 1)).sqrt() * i_n)
    rad = 0.0
    Am = arb_mat(N, N)
    for i in range(N):
        for j in range(N):
            v = G[i, j] + sign * 2 * pole[i] * pole[j] + (beta if i == j else 0)
            Am[i, j] = v.mid()
            rad = max(rad, float(v.rad()))
    vals = smallest_two(Am, job)
    for kb in [int(x) for x in args.blocks.split(",") if x]:
        if kb < N:
            sub = arb_mat([[Am[i, j] for j in range(kb)] for i in range(kb)])
            lb = smallest_two(sub, job)[0]
            job.result(f"leading block k={kb} (deg ≤ {degrees[kb - 1]}): lambda_1 = {mp.nstr(lb, 15)}")
    job.result(f"lambda_1 = {mp.nstr(vals[0], 15)}, lambda_2 = {mp.nstr(vals[1], 8)}; max entry radius {rad:.2e}; "
               f"beta = {beta.str(8)}; nodes {len(nodes)}, modes {N} (deg ≤ {degrees[-1]})")
    job.done()


if __name__ == "__main__":
    main()
