"""Theorem (computer-assisted). For every L in (0, L_max], the lowest eigenvalue of Weil's quadratic form on
L^2(-L, L) is simple and its eigenfunction is even.

Piece II, monotonicity. Enlarging the window adds test functions, and Q itself does not change (a prime power n only
acts through g(log n), which vanishes when log n >= 2L). By min-max, every sector eigenvalue is nonincreasing in L.
So on [L_i, L_{i+1}]:  lam0(L) <= U(L_i)  and  lam1even(L), lam0odd(L) >= lower bounds at L_{i+1}.
Piece I, small L (lipschitz.py). With L <= L_0 = 0.05 < log(2)/2 there are no primes. Gaps of Q_L equal gaps of the
rescaled R_L, and those are Lipschitz in L. So for L <= L_0:
    gap_even(L) >= lam1even_lo(L_0) - U(L_0) - 2 C_even L_0,  gap_odd(L) >= lam0odd_lo(L_0) - U(L_0) - (C_even + C_odd) L_0.
Simple + even <=> lam0 < min(lam1even, lam0odd)   (the full spectrum is the union of the two parity sectors).
L_max is the last grid window of the longest certified prefix of GRID (see ../data/theorem_bounds.json).

    python theorem.py           certify missing windows (cache: ./theorem_bounds.json), then glue
    python theorem.py --check   glue only, from ../data/theorem_bounds.json (no window is recomputed; ~1 min)
"""
import json
import os
import sys
import time
from contextlib import nullcontext
from itertools import takewhile
from multiprocessing import Pool
from flint import arb
import certify
import lipschitz

GRID = ["0.05", "0.1", "0.2", "0.3", "0.35", "0.4", "0.45", "0.5", "0.55", "0.6", "0.64", "0.68", "0.72", "0.76", "0.8",
        "0.84", "0.87", "0.9", "0.93", "0.96", "0.97",  # past 0.8: CERT_PREC=1536 (prime 5 enters; T# = 830)
        "0.98", "0.99", "1.0",  # prime 7 enters (T# = 3610): CERT_BESSEL_PREC=6144 CERT_MATMUL_PREC=320 CERT_VCHOL_PREC=384
        "1.03", "1.06", "1.09", "1.12", "1.15", "1.18", "1.19"]  # sharp T# (certify.SHARP_TS); 8, 9 enter; 11 at 1.1989

if __name__ == "__main__":
    t0 = time.time()
    check = "--check" in sys.argv
    bounds = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "theorem_bounds.json") \
        if check else "theorem_bounds.json"
    try:  # windows already certified by this same code are reused
        rows = {L: tuple(arb(x) for x in r) for L, r in json.load(open(bounds)).items()}
    except FileNotFoundError:
        rows = {}
    grid = list(takewhile(lambda L: L in rows, GRID)) if check else GRID
    with nullcontext() if check else Pool(20) as pool:
        for L in grid:
            if L not in rows and float(L) > float(os.environ.get("UPTO", "inf")):
                raise SystemExit(f"stopped before L={L} (UPTO); rerun with the next precision setting")
            if L not in rows:
                rows[L] = certify.window_bounds(L, pool)
                json.dump({k: [x.str(20, radius=True) for x in r] for k, r in rows.items()},
                          open("theorem_bounds.json", "w"), indent=1)
            lo, U, l1, od = rows[L]
            print(f"L={L:5}  {lo.str(4)} <= lam0 <= {U.str(4)}   lam1even >= {l1.str(4)}   lam0odd >= {od.str(4)}"
                  f"   ({time.time() - t0:.0f}s)", flush=True)

    ok = True
    _, C_e, C_o = lipschitz.constants(grid[0])
    lo, U, l1, od = rows[grid[0]]
    L0 = arb(grid[0])
    ge, go = l1 - U - 2 * C_e * L0, od - U - (C_e + C_o) * L0
    print(f"Piece I  (0, {grid[0]}]:  gap_even >= {ge.str(4)}   gap_odd >= {go.str(4)}   {ge > 0 and go > 0}")
    ok &= ge > 0 and go > 0
    for a, b in zip(grid, grid[1:]):
        U, (l1, od) = rows[a][1], rows[b][2:]
        step = min(l1, od) > U
        print(f"Piece II [{a}, {b}]:  min(lam1even, lam0odd)({b}) / lam0_up({a}) >= {(min(l1, od) / U).lower().str(4)}   {step}")
        ok &= step
    print(f"THEOREM CERTIFIED for all L in (0, {grid[-1]}]" if ok else "NOT certified: refine the grid where a step failed")
