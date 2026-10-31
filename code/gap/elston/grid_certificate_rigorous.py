#!/usr/bin/env python3
"""Certified lower bound for Weil's form on [−a, a] using the window-compressed prime norm.

Derivation (docs/GRID_NORM.md). For real f supported in [−a, a], with F(t) = ∫ f e^{itu} du,
  Q(f) = pole + (1/π) ∫₀^∞ Ψ |F|²,  Ψ = H − P,  H(t) = Re ψ(1/4 + it/2) − log π,
  P(t) = Σ_{log n < 2a} (2Λ(n)/√n) cos(t log n).
Let κ = (sin Ωu / πu)·W(u), where W is the Kaiser window I₀(β_K√(1 − u²/δ²))/I₀(β_K) on [−δ, δ].
Put K = κ̂ and m = 1 − K, and set h = f − f∗κ, so that ĥ = mF and supp h ⊂ [−b, b], b = a + δ.
For μ̄ ≥ sup{⟨g, P g⟩/‖g‖² : supp g ⊂ [−b, b]} (a rigorous Schur bound, scripts/grid_norm.py),
  Q(f) ≥ pole + (1/π) ∫₀^∞ Φ |F|²,  Φ = Ψ + m²(P − μ̄).
For t ≥ T♯ = Ω + β_K/δ we have |K| ≤ ε, so Φ ≥ log(t/2π) − 1/t − μ̄ − (2ε + ε²)(A + μ̄) ≥ β′.
Zhu's reduction (Thm 1.1; audited in docs/AUDIT_ZHU.md) then applies with Φ in place of Ψ:
  Q ≥ (min(λ_min(A_N), β′ − ε_D) − ε_B) ‖f‖²,
where A_N is the Legendre block of R′(f) = pole + (1/π)∫₀^{T♯}(Φ − β′)|F|² + β′‖f‖².

Every input to the bound is enclosed in Arb:
- the Gauss–Legendre nodes and weights;
- Ψ, through acb.digamma;
- m: the main lobe by acb.integral, the sidelobes by the second mean value theorem;
- j_n: the ratio enclosure above the turning point, and upward recurrence below it at adaptive precision;
- the Gram sums, the pole vector and the quadrature-truncation bound (Bernstein ellipse);
- the coupling bounds, and the Cholesky residual (computed exactly with integers).

Usage:
  .venv/bin/python scripts/grid_certificate_rigorous.py --a 1.19 --delta 0.3 --beta-k 64 --tsharp 700 \
      --sector even --nmodes 660 --gl 56 --workers 4 --json data/connes/grid_certificate_a119_even.json
"""

import argparse
import json
import math
import multiprocessing as mproc
import os
import sys
import time

from flint import acb, arb, arb_mat, ctx, fmpz, fmpz_mat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from progress import Job  # noqa: E402
import grid_norm as gn  # noqa: E402

PREC = 400


# ---------------------------------------------------------------- constants

def comb_terms(a):
    out = []
    for n in range(2, int(math.exp(2 * a)) + 1):
        if math.log(n) >= 2 * a:
            continue
        m, p = n, next(q for q in range(2, n + 1) if n % q == 0)
        while m % p == 0:
            m //= p
        if m == 1:
            out.append((n, arb(n).log(), 2 * arb(p).log() / arb(n).sqrt()))
    return out


def mu_upper(a, b, K, iters):
    """Rigorous Schur upper bound on the top of the comb's shift operator on L²[−b, b]."""
    w = gn.weights("zeta", b, comb_a=a)
    _, _, _, vec = gn.top_eigenvalue(w, b, K)
    n_list = [n for n in range(2, int(math.exp(2 * b)) + 1) if any(abs(math.log(n) - l) < 1e-12 for l, _ in w)]
    phi = gn.schur_weight(n_list, b, K, vec, iters=iters)
    ctx.prec = PREC
    return gn.schur_upper(n_list, b, K, phi), n_list


# ---------------------------------------------------------------- Bessel enclosures

def sph_j_enclosure(x, nmax):
    """Enclosures of j_0(x) … j_nmax(x) for real x > 0 (x a ball):
    upward recurrence from j₀ = sinc x and j₁ up to the turning point n_tp = ⌊x⌋, at extra precision
    to absorb the ball-arithmetic wrapping; above it, ratios r_n = j_n/j_{n−1} ∈ (0, 1) by downward
    iteration of r_n = 1/((2n+1)/x − r_{n+1}) from r ∈ [0, 1] (as in docs/AUDIT_ZHU.md)."""
    xf = float(x.mid())
    n_tp = min(int(math.floor(xf)), nmax)
    x_rad = x.rad()
    old = ctx.prec
    ctx.prec = PREC + int(1.25 * xf) + 96
    try:
        # Evaluate at the exact dyadic midpoint (radius 0): the recurrence amplifies an input radius
        # by about e^{0.47x}, independent of the working precision. The node's own radius is added
        # back at the end, using |j_n'(x)| ≤ 1/2 (j_n = ½(−i)^n ∫_{−1}^{1} e^{ixs} P_n(s) ds).
        xx = arb(x.mid())
        j = [None] * (nmax + 1)
        j[0] = xx.sinc()
        if nmax >= 1 and n_tp >= 1:
            j[1] = (j[0] - xx.cos()) / xx
            for n in range(1, n_tp):
                j[n + 1] = (2 * n + 1) / xx * j[n] - j[n - 1]
        if nmax > n_tp:
            # The ratio map r ↦ 1/((2n+1)/x − r) contracts above the turning point, so it needs no
            # extra precision; only the upward recurrence below n_tp does (ball wrapping ~e^{0.47x}).
            ctx.prec = PREC + 64
            xr = arb(x.mid())
            n_start = max(int(2 * xf), nmax) + 60
            r = arb(0.5, 0.5)
            ratios = {}
            for n in range(n_start, n_tp, -1):
                r = 1 / ((2 * n + 1) / xr - r)
                if n <= nmax:
                    ratios[n] = r
            for n in range(n_tp + 1, nmax + 1):
                j[n] = j[n - 1] * ratios[n]
    finally:
        ctx.prec = old
    pad = arb(0, (x_rad / 2).upper()) if x_rad != 0 else arb(0)
    return [arb(v) + pad for v in j]


# ---------------------------------------------------------------- workers

_G = {}


def _init():
    ctx.prec = PREC
    ctx.threads = 1


def _node_value(idx):
    t, wt = _G["nodes"][idx]
    L, degrees, comb, mu, beta = _G["L"], _G["degrees"], _G["comb"], _G["mu"], _G["beta"]
    H = acb(arb(1) / 4, t / 2).digamma().real - _G["log_pi"]
    P = sum((c * (t * l).cos() for _, l, c in comb), arb(0))
    m = _G["mtab"][idx]
    phi = (H - P) + m * m * (P - mu)
    cq = wt * (phi - beta) / arb.pi()
    if t == 0:
        jv = [arb(1)] + [arb(0)] * degrees[-1]
    else:
        jv = sph_j_enclosure(L * t, degrees[-1])
    col = [(2 * L * (2 * n + 1)).sqrt() * jv[n] for n in degrees]
    return cq, col


def _chunk_spans(task):
    """Node spans [s, e) of the chunks one worker processes. Every chunk is contiguous in t (the FLINT
    product's cost depends on the within-chunk spread of j_n ∝ xⁿ, so chunks are never interleaved).
      ("range", lo, hi):        contiguous block of nodes, chunked (the original static partition);
      ("strided", wid, nw):     chunks wid, wid+nw, wid+2nw, … of the whole range, so the expensive
                                low-t chunks are dealt round-robin to all workers instead of one."""
    chunk, n = _G["chunk"], len(_G["nodes"])
    if task[0] == "range":
        lo, hi = task[1], task[2]
        return [(s, min(hi, s + chunk)) for s in range(lo, hi, chunk)]
    wid, nw = task[1], task[2]
    nchunks = (n + chunk - 1) // chunk
    return [(c * chunk, min(n, (c + 1) * chunk)) for c in range(wid, nchunks, nw)]


def _mant_bytes():
    """Fixed byte width of a mantissa record: the midpoint mantissa has at most PREC bits (the radius mantissa
    30), plus a sign; sized from PREC so an overflow is impossible at any precision the run uses."""
    return PREC // 8 + 16


def _work(task):
    N = len(_G["degrees"])
    G = arb_mat(N, N)
    for s, e in _chunk_spans(task):
        cols, scaled = [], []
        for idx in range(s, e):
            cq, col = _node_value(idx)
            cols.append(col)
            scaled.append([cq * v for v in col])
        k = len(cols)
        S = arb_mat([[cols[j][i] for j in range(k)] for i in range(N)])
        D = arb_mat([[scaled[j][i] for j in range(k)] for i in range(N)])
        G += D * S.transpose()
        _G["q"].put(k)
    # Write the upper triangle as fixed-width binary records (mid mantissa, mid exponent, radius mantissa,
    # radius exponent: two's-complement little-endian ints). Exact, so the main process rebuilds the same
    # balls; int.to_bytes raises OverflowError rather than truncating if a value does not fit.
    mb = _mant_bytes()
    tag = f"{task[1]}_{task[2]}" if task[0] == "range" else f"w{task[1]}of{task[2]}"
    path = os.path.join(_G["tmpdir"], f"gram_part_{tag}_{os.getpid()}.bin")
    with open(path, "wb") as fh:
        fh.write(N.to_bytes(4, "little") + mb.to_bytes(4, "little"))
        for i in range(N):
            buf = bytearray()
            for j in range(i, N):
                v = G[i, j]
                mm, me = v.mid().man_exp()
                rm, re_ = v.rad().mid().man_exp()
                buf += int(mm).to_bytes(mb, "little", signed=True) + int(me).to_bytes(8, "little", signed=True)
                buf += int(rm).to_bytes(mb, "little", signed=True) + int(re_).to_bytes(8, "little", signed=True)
            fh.write(buf)
    return path


def _from_pack(t):
    mm, me, rm, re_ = t
    return arb(arb(fmpz(mm)) * arb(2) ** me, arb(fmpz(rm)) * arb(2) ** re_)


def _read_part(path, N):
    """Rebuild one worker's symmetric partial Gram matrix from its byte-record file (exact)."""
    fb = int.from_bytes
    with open(path, "rb") as fh:
        n_file, mb = fb(fh.read(4), "little"), fb(fh.read(4), "little")
        assert n_file == N, (n_file, N)
        rec = 2 * mb + 16
        rows = [[None] * N for _ in range(N)]
        for i in range(N):
            buf = fh.read(rec * (N - i))
            assert len(buf) == rec * (N - i), "short partial Gram file"
            for jj in range(N - i):
                o = jj * rec
                v = _from_pack((fb(buf[o:o + mb], "little", signed=True), fb(buf[o + mb:o + mb + 8], "little", signed=True),
                                fb(buf[o + mb + 8:o + 2 * mb + 8], "little", signed=True), fb(buf[o + 2 * mb + 8:o + rec], "little", signed=True)))
                rows[i][i + jj] = v
                rows[i + jj][i] = v
    return arb_mat(rows)


_M = {}


def _mtab_segment(seg):
    """m(t) for a contiguous run of sorted band nodes: one integral from the lobe start to the first
    node, then cumulative integrals between consecutive nodes (all acb.integral, single-threaded)."""
    ctx.prec = PREC
    ctx.threads = 1
    nodes, what, lo_e, Omega = _M["nodes"], _M["what"], _M["lo_e"], _M["Omega"]
    out = []
    Gcum, prev = arb(0), lo_e
    for _, idx in seg:
        t = nodes[idx][0]
        tau_e = (t.mid() - Omega).mid()
        Gcum += acb.integral(what, acb(prev), acb(tau_e)).real
        prev = tau_e
        m = _M["c_left"] + Gcum + _M["slack"] + kaiser_right_tail(t + Omega, _M["bk"], _M["delta"], _M["i0"])
        mm, me = m.mid().man_exp()
        rm, re_ = m.rad().mid().man_exp()
        out.append((idx, (int(mm), int(me), int(rm), int(re_))))
    return out


# ---------------------------------------------------------------- certification helpers

def cholesky_mid(A):
    """Floating (midpoint) Cholesky factor of a symmetric arb_mat, as nested lists of exact dyadics."""
    N = A.nrows()
    Lm = [[arb(0)] * N for _ in range(N)]
    for j in range(N):
        s = A[j, j].mid() - sum((Lm[j][k] * Lm[j][k] for k in range(j)), arb(0)).mid()
        if not s > 0:
            return None, j
        d = s.sqrt().mid()
        Lm[j][j] = d
        Lj = Lm[j]
        for i in range(j + 1, N):
            Li = Lm[i]
            v = A[i, j].mid() - sum((Li[k] * Lj[k] for k in range(j)), arb(0)).mid()
            Li[j] = (v / d).mid()
    return Lm, None


def cholesky_blocked(A, bs=128):
    """Midpoint Cholesky of a symmetric arb_mat by blocks: bs×bs diagonal factors in Python,
    panel solves (arb_mat.solve) and trailing updates (arb_mat products) in FLINT. Same output as
    cholesky_mid (nested lists of exact dyadics) for the exact-residual certificate."""
    N = A.nrows()
    Lm = [[arb(0)] * N for _ in range(N)]
    for k0 in range(0, N, bs):
        k1 = min(N, k0 + bs)
        m = k1 - k0
        D = arb_mat([[A[k0 + i, k0 + j].mid() for j in range(m)] for i in range(m)])
        if k0 > 0:
            Lk = arb_mat([[Lm[k0 + i][j] for j in range(k0)] for i in range(m)])
            D = (D - Lk * Lk.transpose()).mid()
        Ld, fail = cholesky_mid(D)
        if Ld is None:
            return None, k0 + fail
        for i in range(m):
            for j in range(i + 1):
                Lm[k0 + i][k0 + j] = Ld[i][j]
        if k1 < N:
            B = arb_mat([[A[k1 + i, k0 + j].mid() for j in range(m)] for i in range(N - k1)])
            if k0 > 0:
                Lb = arb_mat([[Lm[k1 + i][j] for j in range(k0)] for i in range(N - k1)])
                B = (B - Lb * Lk.transpose()).mid()
            Lkk = arb_mat([[Ld[i][j] if j <= i else arb(0) for j in range(m)] for i in range(m)])
            Xt = Lkk.solve(B.transpose(), nonstop=True).mid()      # Lkk · Xᵀ = Bᵀ  ⇒  X = B · Lkk^{−T}
            for i in range(N - k1):
                row = Lm[k1 + i]
                for j in range(m):
                    row[k0 + j] = Xt[j, i]
    return Lm, None


def _floor_scaled(x, s):
    """floor(x · 2^s) exactly, for an exact dyadic arb x (its midpoint)."""
    man, ex = x.mid().man_exp()
    man, e = int(man), int(ex) + s
    return man << e if e >= 0 else man >> (-e)  # Python's >> floors, also for negatives


def exact_residual(A, Lm, mu, P):
    """Exact Gram residual (audit fix: no rounded arithmetic before the integer step).
    With L̃ = Lint/2^P, μ′ = floor(μ_lower·2^{2P})/2^{2P} ≤ μ and Aint = floor(A_mid·2^{2P}),
    all exact integers, R = Aint − μ′·2^{2P} I − Lint Lintᵀ is exact. Each floor of A_mid moves an entry
    by < 1 unit, which costs N units per row, so λ_min(A_mid) ≥ μ′ − ‖R‖_∞/2^{2P} − N/2^{2P}.
    Returns (residual bound, μ′)."""
    N = A.nrows()
    Lint = fmpz_mat(N, N, [_floor_scaled(Lm[i][j], P) for i in range(N) for j in range(N)])
    mu_int = _floor_scaled(mu.lower(), 2 * P)
    Aint = fmpz_mat(N, N, [_floor_scaled(A[i, j], 2 * P) - (mu_int if i == j else 0) for i in range(N) for j in range(N)])
    R = Aint - Lint * Lint.transpose()
    worst = 0
    for i in range(N):
        worst = max(worst, sum(abs(int(R[i, j])) for j in range(N)))
    s2 = arb(2) ** (2 * P)
    return (arb(worst) + N) / s2, arb(mu_int) / s2


def coupling_bounds(degrees, L, T, Greal, pole_abs, nmax_tail):
    """ε_B and ε_D as in docs/AUDIT_ZHU.md, with G = sup_{[0,T♯]} |Φ − β′|."""
    n0 = degrees[-1] + 2

    def dfact(n):
        v = arb(1)
        for k in range(1, 2 * n + 2, 2):
            v *= k
        return v

    def d(n):
        return (2 * L * (2 * n + 1)).sqrt() * (L * T) ** n / dfact(n)

    def e(n):
        z = L / 2
        return (2 * L * (2 * n + 1)).sqrt() * z ** n / dfact(n) * (z * z / (2 * (2 * n + 3))).exp()

    def q_of(n, y):
        return (arb(2 * n + 5) / (2 * n + 1)).sqrt() * y ** 2 / ((2 * n + 3) * (2 * n + 5))

    n1 = n0
    while float(q_of(n1, L * T).upper()) >= 0.5:
        n1 += 2
    q_d, q_e = q_of(n1, L * T), q_of(n1, L / 2)
    sum_d = sum((d(n) for n in range(n0, n1, 2)), arb(0)) + d(n1) / (1 - q_d)
    sum_e = sum((e(n) for n in range(n0, n1, 2)), arb(0)) + e(n1) / (1 - q_e)
    c = Greal * T / arb.pi()
    s_inner = [(2 * L * (2 * m + 1)).sqrt() for m in degrees]
    B_inf = 2 * max(pole_abs, key=lambda z: float(z.mid())) * sum_e + c * sum_d * max(s_inner, key=lambda z: float(z.mid()))
    B_one = 2 * sum_e * sum(pole_abs, arb(0)) + c * sum_d * sum(s_inner, arb(0))
    eps_B = (B_inf * B_one).sqrt()
    e0, d0 = e(n0), d(n0)
    eps_D = 2 * e0 * sum_e + c * d0 * sum_d
    return {"n0": n0, "d_n0": d0, "eps_B": eps_B, "eps_D": eps_D}


def kaiser_right_tail(s0, bk, delta, i0, K=40):
    """(1/2π)∫_{s0}^∞ Ŵ(s) ds for s0 > β_K/δ, where Ŵ is the Kaiser window's Fourier transform.
    With r = √(δ²s² − β²) the integral is (1/(π I₀)) ∫_{r0}^∞ sin(r) g(r) dr, g = (r² + β²)^{-1/2}.
    Repeated integration by parts gives ∫_{r0}^∞ e^{ir} g = e^{i r0} Σ_{k<K} i^{k+1} g^{(k)}(r0) + i^K ∫ e^{ir} g^{(K)},
    and |g^{(K)}(r)| ≤ K!/r^{K+1} (g(r) = ∫₀^∞ J₀(βs) e^{−rs} ds, |J₀| ≤ 1), so the remainder is at
    most (K−1)!/r0^K. The Taylor coefficients of g come from the power recurrence for (p0 + p1x + x²)^{−1/2}."""
    r0 = ((delta * s0) ** 2 - bk ** 2).sqrt()
    p0, p1 = r0 * r0 + bk * bk, 2 * r0
    alpha = arb(-1) / 2
    f = [p0 ** alpha]
    for n in range(1, K):
        acc = ((alpha + 1) * 1 - n) * p1 * f[n - 1]
        if n >= 2:
            acc += ((alpha + 1) * 2 - n) * f[n - 2]
        f.append(acc / (n * p0))
    S = acb(0)
    fac = arb(1)
    ipow = acb(0, 1)
    for k in range(K):
        if k > 0:
            fac *= k
        S += ipow * (fac * f[k])
        ipow *= acb(0, 1)
    I = acb(r0.cos(), r0.sin()) * S
    rem = fac * K / K / r0 ** K  # (K−1)!/r0^K, fac = (K−1)!
    val = I.imag + arb(0, rem.upper())
    return val / (arb.pi() * i0)


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--a", default="1.19")
    ap.add_argument("--delta", default="0.3")
    ap.add_argument("--beta-k", default="64")
    ap.add_argument("--tsharp", default="700")
    ap.add_argument("--sector", choices=["even", "odd"], default="even")
    ap.add_argument("--nmodes", type=int, default=660)
    ap.add_argument("--gl", type=int, default=56)
    ap.add_argument("--panels-per-unit", type=int, default=4, help="panel width 1/this (ρ = 2+√5 needs width 1/4)")
    ap.add_argument("--schur-k", type=int, default=4000)
    ap.add_argument("--schur-iters", type=int, default=1500)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--chunk", type=int, default=512)
    ap.add_argument("--json")
    ap.add_argument("--plan", action="store_true", help="compute every error budget, skip the Gram matrix")
    ap.add_argument("--prec", type=int, default=400, help="working precision in bits (Bessel recurrences add more)")
    ap.add_argument("--cholesky", choices=["blocked", "python"], default="blocked")
    ap.add_argument("--block", type=int, default=128, help="block size for the blocked Cholesky")
    ap.add_argument("--tmpdir", default=None, help="directory for the workers' partial Gram files (default: logs/)")
    ap.add_argument("--partition", choices=["strided", "contiguous"], default="strided",
                    help="how chunks are dealt to workers: strided (round-robin, no straggler) or contiguous (the original blocks)")
    ap.add_argument("--deflate", action="store_true", help="GenusCore: certify λ₂ (even sector, simple-even) via rank-one deflation")
    ap.add_argument("--lambda2", action="store_true", help="also compute λ₂ (diagnostic gap; about 40 min at N = 3200)")
    args = ap.parse_args()
    global PREC
    PREC = args.prec
    ctx.prec = PREC
    ctx.threads = 1  # acb.integral calls back into Python: FLINT threads would segfault
    t_start = time.time()
    job = Job(f"grid certificate rigorous {args.sector} a={args.a} T{args.tsharp} N{args.nmodes}", args=vars(args))
    out = {"args": vars(args)}

    a, delta, bk, T = arb(args.a), arb(args.delta), arb(args.beta_k), arb(args.tsharp)
    L = a
    from fractions import Fraction
    b_exact = Fraction(args.a) + Fraction(args.delta)
    b = float(b_exact) + 1e-9  # Schur window, rounded up: μ is monotone in the window
    assert Fraction(b) >= b_exact
    comb = comb_terms(float(args.a))
    A_sup = sum((c for _, _, c in comb), arb(0))
    log_pi = arb.pi().log()

    # 1. μ̄ (Schur)
    mu, n_list = mu_upper(float(args.a), b, args.schur_k, args.schur_iters)
    mu = arb(mu)
    job.log(f"comb n = {[n for n, _, _ in comb]}; A = {A_sup.str(12)}; window b = {b}; Schur μ̄ = {mu.str(25)} (shifts {n_list})")
    out.update({"comb": [n for n, _, _ in comb], "A": A_sup.str(20), "b": b, "mu_upper": mu.str(20)})

    # 2. Kaiser kernel and β′
    half = bk / delta
    Omega = T - half
    i0 = acb(bk).bessel_i(0).real
    eps_side = 2 / (arb.pi() * bk * i0)
    assert float((Omega * delta).lower()) >= 1 and float((Omega - half).lower()) > 0
    wbound = 2 * eps_side + eps_side ** 2
    beta_bound = (T / (2 * arb.pi())).log() - 1 / T - mu - wbound * (A_sup + mu)
    beta = (((beta_bound.lower() * arb(2) ** 60).floor() - 1) / arb(2) ** 60).mid()  # exact dyadic, strictly below
    if not (beta < beta_bound and beta > 0):
        raise SystemExit(f"β′ check failed: beta = {beta}, bound = {beta_bound}, mu = {mu}, eps = {eps_side}")
    job.log(f"Ω = {Omega.str(10)}, transition [{(Omega - half).str(8)}, {T.str(8)}]; sidelobe bound ε = {eps_side.str(4)}; β′ = {beta.str(15)} (bound {beta_bound.str(12)})")
    out.update({"Omega": Omega.str(20), "eps_side": eps_side.str(6), "beta_prime": beta.str(25)})

    # 3. nodes (exact GL balls)
    gl = [arb.legendre_p_root(args.gl, i, weight=True) for i in range(args.gl)]
    panels = int(args.tsharp) * args.panels_per_unit
    h = T / panels
    nodes = []
    for p in range(panels):
        c = h * p + h / 2
        for x, wt in gl:
            nodes.append((c + x * h / 2, wt * h / 2))
    job.log(f"{len(nodes)} nodes ({panels} panels × {args.gl}-point Gauss)")

    # 4. m(t): passband/sidelobe bound, main lobe by rigorous integration
    pref = delta / (arb.pi() * i0)  # Ŵ(s)/(2π) = pref · 0F1(;3/2;(β² − δ²s²)/4)

    def what(z, analytic):
        return acb(pref) * ((acb(bk) ** 2 - (acb(delta) * z) ** 2) / 4).hypgeom_0f1(acb(1.5))

    w_peak = pref * (bk.sinh() / bk)
    lo_band = Omega - half
    mtab = [None] * len(nodes)
    band = []
    for idx, (t, _) in enumerate(nodes):
        if t < lo_band:
            mtab[idx] = arb(0, 2 * eps_side)
        elif t > lo_band:
            band.append((float(t.mid()), idx))
        else:
            raise ValueError("node overlaps the transition start; shift T♯")
    band.sort()
    # Exact dyadic integration endpoints. The true endpoints (−β_K/δ and t − Ω, both balls) differ
    # from them by less than 1e-100, which moves the integral by less than 1e-100·w_peak; the split
    # point between "tail" and "main lobe" is free, and a sliver of main lobe is covered by that radius.
    slack = arb(0, (arb("1e-100") * w_peak).upper())
    # m(t) = c_left + M(t) + R(t) in the band, where c_left = (1 − main)/2 by symmetry of Ŵ,
    # M(t) = (1/2π)∫_{−β_K/δ}^{t−Ω} Ŵ, and R(t) = (1/2π)∫_{t+Ω}^∞ Ŵ (kaiser_right_tail).
    lo_e, hi_e = (-half).mid(), half.mid()
    main_total = acb.integral(what, acb(lo_e), acb(hi_e)).real
    c_left = (1 - main_total) / 2 + slack
    job.log(f"Kaiser main lobe (1/2π)∫Ŵ = {main_total.str(20)}; left tail c = {c_left.str(6)}")
    nw_m = max(1, min(args.workers, len(band) // 2000 or 1))
    segs = [band[len(band) * i // nw_m: len(band) * (i + 1) // nw_m] for i in range(nw_m)]
    _M.update({"nodes": nodes, "what": what, "lo_e": lo_e, "Omega": Omega, "c_left": c_left, "slack": slack,
               "bk": bk, "delta": delta, "i0": i0})
    t_m = time.time()
    if nw_m == 1:
        results = [_mtab_segment(segs[0])]
    else:
        with mproc.get_context("fork").Pool(nw_m, initializer=_init) as pool:
            results = pool.map(_mtab_segment, segs)
    ctx.prec = PREC
    for res_seg in results:
        for idx, tpl in res_seg:
            mtab[idx] = _from_pack(tpl)
    job.log(f"m-table: {len(band)} band nodes in {nw_m} segments, {time.time() - t_m:.0f}s")
    job.log(f"m at the last band node = {mtab[band[-1][1]].str(12)}")

    first = 0 if args.sector == "even" else 1
    degrees = [first + 2 * i for i in range(args.nmodes)]
    sign = 1 if args.sector == "even" else -1
    N = len(degrees)
    # 5. pole vector
    pole = []
    for n in degrees:
        sg = (-1) ** (n // 2) if args.sector == "even" else (-1) ** ((n - 1) // 2)
        xh = L / 2
        i_n = (arb.pi() / (2 * xh)).sqrt() * acb(xh).bessel_i(acb(n + arb(1) / 2)).real
        pole.append(sg * (2 * L * (2 * n + 1)).sqrt() * i_n)
    # 7. quadrature truncation (Bernstein ellipse ρ = 2+√5 on panels of width 1/4: semi-minor 1/4)
    rho = 2 + arb(5).sqrt()
    b_e = h * (rho - 1 / rho) / 4
    assert float(b_e.upper()) <= 0.25 + 1e-12
    re_w = arb(1) / 4 - b_e / 2
    v_re = 1 + re_w
    v_abs = 1 + arb(1) / 4 + b_e / 2 + (T + h * (rho + 1 / rho) / 4) / 2
    psi_abs = 1 / re_w + v_abs.log() + arb.pi() / 2 + 1 / (2 * v_re) + 1 / (12 * v_re ** 2)
    comb_abs = sum((c * (b_e * l).cosh() for _, l, c in comb), arb(0))
    m_abs = 1 + (delta * b_e).exp() * (2 / arb.pi()) * (1 + (Omega * delta).log())
    G_strip = psi_abs + log_pi + comb_abs + m_abs ** 2 * (comb_abs + mu) + beta
    S2 = 2 * L * (2 * degrees[-1] + 1) * (2 * L * b_e).exp()
    Mb = G_strip * S2 / arb.pi()
    eps_Q = T * Mb * 4 * rho ** (-2 * args.gl) / (1 - 1 / rho)
    # Sharper spectral debit (audit, Cauchy–Schwarz): entry errors obey |E_ij| ≤ C u_i u_j, where
    # u_i = √(2L(2n_i+1)) and C = ε_Q/u_max², so ‖E‖₂ ≤ C Σ u_i² ≤ N ε_Q.
    u2_sum = sum((2 * L * (2 * n + 1) for n in degrees), arb(0))
    c_err = eps_Q / (2 * L * (2 * degrees[-1] + 1)) * u2_sum
    job.log(f"quadrature: |Φ−β′| ≤ {G_strip.str(6)} on the strip, M ≤ {Mb.str(6)}, ε_Q ≤ {eps_Q.str(4)}, c_err = N ε_Q ≤ {c_err.str(4)}")

    # 8. coupling (real axis): |H| ≤ 8.1 on [0, T♯] (ψ(1/4) ≤ Re ψ ≤ log|z| + 2 + 4/3), m ≤ 1 + 2ε
    Re_psi_max = (arb(1) / 16 + T * T / 4).sqrt().log() + 2 + arb(4) / 3
    H_abs = max(abs(acb(arb(1) / 4).digamma().real - log_pi), Re_psi_max - log_pi, key=lambda z: float(z.upper()))
    G_real = H_abs + A_sup + (1 + 2 * eps_side) ** 2 * (A_sup + mu) + beta
    coup = coupling_bounds(degrees, L, T, G_real, [abs(p) for p in pole], degrees[-1])
    job.log(f"coupling: G_real = {G_real.str(6)}, n0 = {coup['n0']}, d_n0 = {coup['d_n0'].str(4)}, ε_B = {coup['eps_B'].str(4)}, ε_D = {coup['eps_D'].str(4)}")

    if args.plan:
        total = c_err + coup["eps_B"]
        job.result(f"PLAN: c_err ≤ {c_err.str(4)}, ε_B ≤ {coup['eps_B'].str(4)}, ε_D ≤ {coup['eps_D'].str(4)}; total budget {total.str(4)}")
        job.done()
        return
    # 8. Gram matrix (parallel)
    mgr = mproc.Manager()
    q = mgr.Queue()
    tmpdir = args.tmpdir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs")
    os.makedirs(tmpdir, exist_ok=True)
    cfg = {"nodes": nodes, "L": L, "degrees": degrees, "comb": comb, "mu": mu, "beta": beta,
           "log_pi": log_pi, "mtab": mtab, "chunk": args.chunk, "q": q, "tmpdir": tmpdir}
    nw = args.workers
    if args.partition == "contiguous":
        tasks = [("range", len(nodes) * i // nw, len(nodes) * (i + 1) // nw) for i in range(nw)]
    else:
        tasks = [("strided", i, nw) for i in range(nw)]
    job.total = len(nodes)
    t0 = time.time()
    ctxp = mproc.get_context("fork")
    _G.update(cfg)  # inherited by the forked workers
    with ctxp.Pool(nw, initializer=_init) as pool:
        res = pool.map_async(_work, tasks)
        done = 0
        while not res.ready() or not q.empty():
            try:
                k = q.get(timeout=5)
                done += k
                job.step(f"{done}/{len(nodes)} nodes", n=k)
            except Exception:
                pass
        parts = res.get()
    ctx.prec = PREC
    ctx.threads = args.workers  # matrix work only from here on
    t_c = time.time()
    G = arb_mat(N, N)
    for path in parts:          # one exact rebuild and one matrix add per worker (entrywise, same order as before)
        G += _read_part(path, N)
        os.remove(path)
    job.log(f"Gram assembled in {time.time() - t0:.0f}s (collection {time.time() - t_c:.0f}s, partition {args.partition})")

    A = arb_mat(N, N)
    max_rad = 0.0
    for i in range(N):
        for j in range(N):
            v = G[i, j] + sign * 2 * pole[i] * pole[j] + (beta if i == j else 0)
            A[i, j] = v
            max_rad = max(max_rad, float(v.rad()))
    job.log(f"block assembled; max entry radius {max_rad:.3e}")

    # 9. eigenvalue (numerical) and Cholesky certificate (rigorous)
    from grid_certificate_flint import smallest_two
    Am = arb_mat(N, N)
    for i in range(N):
        for j in range(N):
            Am[i, j] = A[i, j].mid()
    if args.deflate:
        # GenusCore addition (simple-even certificates): certify λ₂ instead of λ₁. Rank-one interlacing:
        # A + v vᵀ − μI ⪰ 0 for any v  ⇒  A has at most one eigenvalue below μ, so λ₂(A) ≥ μ; the block lemma then gives
        # λ₂(R′) ≥ min(λ₂(A), β′ − ε_D) − ε_B, and λ₂(Q) ≥ λ₂(R′) by min-max since Q ≥ R′ as forms.
        # v: unit eigenvector of λ₁(A_mid) by inverse iteration, parsed from 20-digit decimal strings as balls; A + v vᵀ is formed in balls, so
        # its rounding enters max_rad below exactly like the Gram radii.
        v = arb_mat([[arb(1) / (1 + i)] for i in range(N)])
        for _ in range(12):
            v = Am.solve(v, nonstop=True)
            nrm = sum((v[i, 0].mid() ** 2 for i in range(N)), arb(0)).sqrt()
            v = arb_mat([[arb((v[i, 0] / nrm).mid().str(20, radius=False))] for i in range(N)])
        A = arb_mat([[A[i, j] + v[i, 0] * v[j, 0] for j in range(N)] for i in range(N)])
        max_rad = max(float(A[i, j].rad()) for i in range(N) for j in range(N))
        Am = arb_mat([[A[i, j].mid() for j in range(N)] for i in range(N)])
        job.log(f"deflated by v vᵀ (‖v‖ ≈ 1); max entry radius now {max_rad:.3e}; certifying λ₂")
    lam = smallest_two(Am, job, count=2 if args.lambda2 else 1)
    if args.lambda2:
        job.result(f"lambda_1(A_mid) = {lam[0]}, lambda_2 = {lam[1]} (inverse iteration: smallest |λ|)")
    else:
        job.result(f"lambda_1(A_mid) = {lam[0]} (inverse iteration: smallest |λ|; λ₂ skipped, pass --lambda2 for the gap)")
    import mpmath as mp
    if N <= 160:
        mp.mp.prec = PREC
        Mm = mp.matrix(N, N)
        for i in range(N):
            for j in range(N):
                Mm[i, j] = mp.mpf(Am[i, j].str(130, radius=False))
        ev = sorted(mp.eigsy(Mm, eigvals_only=True))
        job.result("full spectrum check (small N): lowest = " + ", ".join(mp.nstr(e, 6) for e in ev[:4]))
    mu_s = arb(str(mp.nstr(lam[0] * (1 - mp.mpf(10) ** -6), 20)))
    Ash = arb_mat(N, N)
    for i in range(N):
        for j in range(N):
            Ash[i, j] = Am[i, j] - (mu_s if i == j else 0)
    t1 = time.time()
    Lm, fail = cholesky_blocked(Ash, args.block) if args.cholesky == "blocked" else cholesky_mid(Ash)
    if Lm is None:
        job.warn(f"Cholesky failed at pivot {fail}")
        lam_cert = None
    else:
        res, mu_prime = exact_residual(Am, Lm, mu_s, P=PREC + 200)
        lam_A = mu_prime - res - N * arb(max_rad) - c_err
        lam_cert = min(lam_A, beta - coup["eps_D"], key=lambda z: float(z.lower())) - coup["eps_B"]
        job.log(f"Cholesky at μ_s = {mu_s.str(12)} in {time.time() - t1:.0f}s; exact residual ‖R‖_∞ ≤ {res.str(4)}")
        job.result(f"CERTIFIED{' λ₂' if args.deflate else ''} (sector {args.sector}, a = {args.a}): Q(f) ≥ {lam_cert.lower().str(20)} ‖f‖² (lower endpoint; truncate, do not round, when quoting)   "
                   f"[λ_A ≥ {lam_A.lower().str(10)}, β′ − ε_D = {(beta - coup['eps_D']).str(8)}, ε_B = {coup['eps_B'].str(3)}]")
    out.update({"deflate": bool(args.deflate), "sector": args.sector, "nmodes": N, "max_degree": degrees[-1], "nodes": len(nodes),
                "lambda_1_mid": str(lam[0]), "lambda_2_mid": str(lam[1]) if args.lambda2 else None, "mu_shift": mu_s.str(20),
                "max_entry_radius": max_rad, "eps_Q": eps_Q.str(6), "c_err": c_err.str(6),
                "G_strip": G_strip.str(8), "G_real": G_real.str(8),
                "eps_B": coup["eps_B"].str(6), "eps_D": coup["eps_D"].str(6), "n0": coup["n0"],
                "certified_lower_bound": lam_cert.lower().str(20) if lam_cert is not None else None,
                "seconds": round(time.time() - t_start)})
    if args.json:
        os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=2)
    job.done()


if __name__ == "__main__":
    main()
