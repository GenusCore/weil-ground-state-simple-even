#!/usr/bin/env python3
"""Export the real head-matrix certificates for L = 0.05 and 0.1 to Lean (untrusted producer).

Needs the Weil code (`weil-mechanism/rerun`, set WEIL_RERUN) with its head-matrix cache and python-flint.
Reads the cached Arb balls of the head matrices A (even, odd), rounds every ball OUTWARD to a box
|A_ij - C_ij / 2^64| <= rad / 2^64 with integer C_ij, picks the deflation vector v and the lower bounds mu,
computes an approximate Cholesky factor R (also untrusted, rounded to 2^-32), runs the exact integer
check of checker.py, and writes SpectralCertificate/Data/W*.lean plus certificates.json.

Nothing computed here is trusted: Lean re-checks every entry of the residual by kernel evaluation.
What IS trusted is the input (that the Arb balls enclose the true head matrix) and the constants
beta, eD, eB (Lemmas (E), (R), (S) of paper/NOTE.md), which are exported as rational bounds.
"""
import hashlib
import json
import math
import os
import sys
from fractions import Fraction as Q
from pathlib import Path

import checker

HERE = Path(__file__).resolve().parent
RERUN = Path(os.environ.get("WEIL_RERUN", HERE.parents[3] / "weil-mechanism" / "rerun"))
BOUNDS = HERE.parents[1] / "data" / "theorem_bounds.json"
WINDOWS = ("0.05", "0.1")


def exact(x):  # exact rational value of an exact arb (midpoint or radius)
    m, e = x.man_exp()
    return Q(int(m)) * Q(2) ** int(e)


def ball(x):
    return exact(x.mid()), exact(x.rad())


def up(q, digits=3):  # round a positive rational up to `digits` significant decimals
    k = digits - 1 - math.floor(math.log10(q))
    return Q(math.ceil(q * Q(10) ** k), 1) / Q(10) ** k if k >= 0 else Q(math.ceil(q / Q(10) ** -k)) * Q(10) ** -k


def down(q, places=30):
    return Q(math.floor(q * 10 ** places), 10 ** places)


def stored(Ldec):
    row = checker.read_json(BOUNDS)[Ldec]
    return [checker.ball(s) for s in row]  # (lo, hi) exact


def center(A):
    """Integer center C (scale 2^64) and one radius covering every ball, rounded outward."""
    n, C, rad = len(A), [], Q(0)
    for i in range(n):
        row = []
        for j in range(i + 1):
            mid, r = ball(A[i][j])
            mid2, r2 = ball(A[j][i])
            assert (mid, r) == (mid2, r2), "head matrix is not symmetric"
            c = round(mid * checker.SCALE)
            row.append(c)
            rad = max(rad, abs(mid - Q(c, checker.SCALE)) + r)
        C.append(row)
    return C, math.ceil(rad * checker.SCALE)


def factor(C, v, M, n):
    """Approximate Cholesky of C + v v^T - M I at 256 bits, rounded to integers at scale 2^32. Untrusted."""
    from flint import arb, ctx
    with ctx.workprec(256):
        S = [[(arb(C[max(i, j)][min(i, j)] + (v[i] * v[j] if v else 0) - (M if i == j else 0)) / checker.SCALE).mid()
              for j in range(n)] for i in range(n)]
        R = [[arb(0)] * n for _ in range(n)]
        for j in range(n):
            s = (S[j][j] - sum(R[j][k] ** 2 for k in range(j))).mid()
            assert s > 0, "shift too large: not positive definite"
            R[j][j] = s.sqrt().mid()
            for i in range(j + 1, n):
                R[i][j] = ((S[i][j] - sum(R[i][k] * R[j][k] for k in range(j))) / R[j][j]).mid()
        return [[round(exact(R[i][k]) * checker.HALF) for k in range(i + 1)] for i in range(n)]


def lowest(C, v, n):
    import certify
    from flint import arb
    A = [[arb(C[max(i, j)][min(i, j)] + (v[i] * v[j] if v else 0)) / checker.SCALE for j in range(n)] for i in range(n)]
    est, vec = certify.lowest(A)
    return exact(est), vec


def certificate(name, C, rad, v, mu, n):
    est, vec = lowest(C, v, n)
    assert mu < est, f"{name}: requested mu {float(mu)} is not below the estimate {float(est)}"
    mu_num = math.ceil(mu * checker.SCALE)
    M = math.floor((mu + (est - mu) / 2) * checker.SCALE)  # shift halfway between mu and the estimate
    R = factor(C, v, M, n)
    # Drop factor entries below 2^(t-32) for the largest t that still passes: the kernel then skips them.
    for t in range(24, -1, -1):
        Rt = [[x if abs(x) >> t else 0 for x in r[:-1]] + r[-1:] for r in R]
        cert = dict(n=n, C=C, Rc=checker.compress(Rt), v=v, rad=rad, M=M, K=0, mu=mu_num)
        cert["K"] = checker.residual_bound(cert)
        if mu_num + n * cert["K"] <= M:
            break
    checker.check(cert)
    kept = sum(len(vals) for _, vals in cert["Rc"])
    print(f"  {name}: mu = {mu_num / checker.SCALE:.6f}  est = {float(est):.6f}  drop |R| < 2^{t - 32}, "
          f"{kept} of {n * (n + 1) // 2} factor entries stored  n K / 2^64 = {n * cert['K'] / checker.SCALE:.3e}  "
          f"slack = {(M - mu_num - n * cert['K']) / checker.SCALE:.4f}")
    return cert, vec


def export(Ldec):
    sys.path.insert(0, str(RERUN))
    import certify
    W, Ae, eDe, eBe = certify.head(Ldec, 0, None)
    Wo, Ao, eDo, eBo = certify.head(Ldec, 1, None)
    n = len(Ae)
    row = stored(Ldec)
    beta = down(exact(W.beta.lower()))
    eD = {0: up(exact(eDe.upper())), 1: up(exact(eDo.upper()))}
    eB = {0: up(exact(eBe.upper())), 1: up(exact(eBo.upper()))}

    # mu: smallest multiple of 2^-64 (eB rounded up) with min(mu, beta - eD) - eB >= the stored upper endpoint, so Lean's bound implies the row.
    def mu_for(col, parity):
        return row[col][1] + eB[parity]

    Ce, rade = center(Ae)
    Co, rado = center(Ao)
    print(f"L = {Ldec}: n = {n}, rad even = {rade / checker.SCALE:.2e}, rad odd = {rado / checker.SCALE:.2e}")
    even, vec = certificate("even", Ce, rade, [], mu_for(0, 0), n)
    v = [round(exact(x.mid()) * checker.HALF) for x in vec]
    defl, _ = certificate("deflated", Ce, rade, v, mu_for(2, 0), n)
    odd, _ = certificate("odd", Co, rado, [], mu_for(3, 1), n)
    sources = {f"head_{Ldec}_{p}_{n}.txt": hashlib.sha256((RERUN / "cache" / f"head_{Ldec}_{p}_{n}.txt").read_bytes()).hexdigest()
               for p in (0, 1)}
    consts = dict(beta_lo=beta, eD_even=eD[0], eB_even=eB[0], eD_odd=eD[1], eB_odd=eB[1])
    return dict(even=even, deflated=defl, odd=odd), consts, sources


def points():
    """The exact dyadic window points L~ = arb(L).upper() at 768 bits, as certify.Window uses them."""
    sys.path.insert(0, str(RERUN))
    import certify  # noqa: F401  (sets the working precision)
    from flint import arb
    return {L: str(exact(arb(L).upper())) for L in WINDOWS}


def main():
    manifest = {"scale": "entries of C, rad, M, K, mu are integers / 2^64; R and v are integers / 2^32",
                "windows": {}, "points": points()}
    (HERE / "SpectralCertificate" / "Data").mkdir(exist_ok=True)
    for Ldec in WINDOWS:
        certs, consts, sources = export(Ldec)
        tag = checker.tag(Ldec)
        (HERE / "SpectralCertificate" / "Data" / f"{tag}.lean").write_text(checker.emit_data(tag, certs, consts))
        manifest["windows"][Ldec] = {"head_cache_sha256": sources,
                                     "constants": {k: str(v) for k, v in consts.items()},
                                     "mu": {k: str(Q(c["mu"], checker.SCALE)) for k, c in certs.items()}}
    (HERE / "certificates.json").write_text(json.dumps(manifest, indent=1) + "\n")


if __name__ == "__main__":
    main()
