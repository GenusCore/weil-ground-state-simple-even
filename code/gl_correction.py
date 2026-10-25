"""Repair of the Gauss-Legendre off-by-one (found in review, 2026-10-23), without re-assembling any head matrix.

The certified runs used error <= (64/15) M rho^(-2n)/(rho^2-1) for n-point rules. The correct n-point bound is
rho^(-2n+2), i.e. rho^2 times larger. So:

  Lower bounds (lam0, lam1^even, lam0^odd). Each comes from a head matrix whose entries carried radius q_old. The true
  entries lie within rho^2 q_old, so the true matrix differs from the certified ball by an extra entrywise <=
  (rho^2-1) q_old, i.e. spectral norm <= N (rho^2-1) q_old =: delta. By Weyl every certified head lower bound drops by at
  most delta, and the final bounds min(mu, beta-eD) - eB drop by at most delta.

  Upper bound on lam0. Any trial function gives one. We recompute it from scratch with the FIXED rayleigh_upper, using
  a fresh trial vector: the even ground state of a (non-rigorous) Galerkin matrix.

Then Piece I (fresh Lipschitz constants) and all Piece II steps are re-checked on the corrected numbers.
Writes theorem_bounds_corrected.json.
"""
import json
from math import ceil, e, exp, pi as fpi
from multiprocessing import Pool
from flint import arb, ctx
import certify as c
import convergence as cv
import lipschitz


def delta(Ldec):
    out = []
    for parity in (0, 1):
        W = c.Window(Ldec, parity)
        L, N, TS, beta = W.L, W.N, W.TS, W.beta
        B, pi = c.B, c.pi
        rho = 2 * B + (4 * B * B + 1).sqrt()
        re_w, abs_w_max = arb(1) / 4 - B / 2, arb(1) / 4 + B / 2 + (arb(TS) + 1) / 2
        psi_e = max(abs_w_max.log(), -re_w.log()) + pi / 2 + 1 / (2 * re_w) + 1 / (12 * re_w * re_w)
        comb_e = sum(2 * v / arb(n).sqrt() * (B * c.log(n)).cosh() for n, v in W.lam.items())
        M = (psi_e + pi.log() + comb_e + beta) * 2 * L * (2 * L * B).exp() / pi
        q_old = TS * arb(32) / 15 * M * rho ** (-2 * c.NG) / (rho ** 2 - 1)
        out.append((N * (rho ** 2 - 1) * q_old).upper())
    return out  # [delta_even, delta_odd]


def new_upper(Ldec):
    L = float(Ldec)
    D = min(2 * ceil(e * L * 2 * fpi * exp(2 * L) / 8) + 60, 298)
    Q = cv.galerkin(L, D)
    _, (v0, _) = cv.eig_low(Q, list(range(0, D, 2)), 2)
    v = [x.mid() for x in v0][:150]
    return c.rayleigh_upper(c.Window(Ldec, 0), v).str(60, radius=True)  # arb is not picklable


if __name__ == "__main__":
    ctx.prec = 768
    rows = {L: [arb(x) for x in r] for L, r in json.load(open("theorem_bounds.json")).items()}
    grid = sorted(rows, key=float)
    with Pool(12) as pool:
        U_new = {L: arb(s).upper() for L, s in zip(grid, pool.map(new_upper, grid))}
    corr = {}
    print(" L      delta_even  delta_odd   lam0_lo'        U_old        U_new        lam1even_lo'   lam0odd_lo'")
    for L in grid:
        lo, U_old, l1, od = rows[L]
        de, do = delta(L)
        lo2, l12, od2 = lo - de, l1 - de, od - do
        U2 = U_new[L]
        assert lo2 <= U2, L
        corr[L] = (lo2, U2, l12, od2)
        print(f"{L:5} {de.str(2):>10} {do.str(2):>10}  {lo2.str(4):>14} {U_old.str(4):>12} {U2.str(4):>12} {l12.str(4):>14} {od2.str(4):>14}")

    ok = True
    _, C_e, C_o = lipschitz.constants(grid[0])
    lo, U, l1, od = corr[grid[0]]
    L0 = arb(grid[0])
    ge, go = l1 - U - 2 * C_e * L0, od - U - (C_e + C_o) * L0
    print(f"Piece I (0, {grid[0]}]: C_even <= {C_e.str(6)}, C_odd <= {C_o.str(6)}; gap_even >= {ge.str(4)}, gap_odd >= {go.str(4)}")
    ok &= ge > 0 and go > 0
    for a, b in zip(grid, grid[1:]):
        r = min(corr[b][2], corr[b][3]) / corr[a][1]
        print(f"Piece II [{a}, {b}]: ratio >= {r.lower().str(4)}  {r > 1}")
        ok &= r > 1
    json.dump({L: [x.str(20, radius=True) for x in r] for L, r in corr.items()},
              open("theorem_bounds_corrected.json", "w"), indent=1)
    print("THEOREM CERTIFIED (corrected) for all L in (0, 1.0]" if ok else "NOT certified after correction")
