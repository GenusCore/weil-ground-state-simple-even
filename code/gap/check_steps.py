"""Piece III gluing: for consecutive windows L_i < L_{i+1}, check min(l2_even, l1_odd)(L_{i+1}) > U(L_i) (paper, Cor. 3.2).
Inputs (directory argument, default ./out): upper*.jsonl (fr_sweep, field lam0_upper, rigorous outward-rounded upper ends; every
run is a valid upper bound, so the smallest one per L is used),
g_<L>_odd.json and g_<L>_l2.json (Elston certifier; field certified_lower_bound, an Arb ball whose lower end is used).
U(1.19): the smaller of Piece II's certified 1.312e-47 (data/theorem_bounds.json) and the K = 640 Galerkin bound in upper640.jsonl
(8.712e-48; the first step needs it).
All comparisons are done in Arb: a step passes only if the lower end of the left side exceeds the upper end of U.
The last line states the range certified by the contiguous passing prefix of steps from 1.19, together with Pieces I-II
(checked separately by code/theorem.py --check on data/theorem_bounds.json).
Usage: python check_steps.py [dir]   (default: data/piece3 of this repository)"""
import sys, json, os
from flint import arb

d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "piece3")
U = {"1.19": arb("1.312e-47")}
for fn in sorted(f for f in os.listdir(d) if f.startswith("upper") and f.endswith(".jsonl")):
    for line in open(os.path.join(d, fn)):
        r = json.loads(line)
        if r.get("par", 0) == 0 and (r["L"] not in U or arb(r["lam0_upper"]) < U[r["L"]]): U[r["L"]] = arb(r["lam0_upper"])
grid = ["1.19", "1.2205", "1.25", "1.2785", "1.3055", "1.3315", "1.3565", "1.3805", "1.4035", "1.4255", "1.447", "1.468", "1.4885", "1.5"]

def lower(L, kind):
    p = os.path.join(d, f"g_{L}_{kind}.json")
    if not os.path.exists(p): return None
    v = json.load(open(p)).get("certified_lower_bound")
    return arb(v).lower() if v else None

ok_all, worst, end, extra = True, None, grid[0], []
print(f"{'step':>17}  {'U(L_i)':>12}  {'l2_even(L_i+1)':>15}  {'l1_odd(L_i+1)':>14}  ratio   pass")
for a, b in zip(grid, grid[1:]):
    u, l2, lo = U.get(a), lower(b, "l2"), lower(b, "odd")
    if None in (u, l2, lo):
        print(f"[{a}, {b}]  missing: {'U ' if u is None else ''}{'l2 ' if l2 is None else ''}{'odd' if lo is None else ''}"); ok_all = False; continue
    m = l2 if l2 < lo else lo
    passed = bool(m > u.upper())
    ratio = m / u.upper()
    if passed and ok_all: end = b; worst = ratio if worst is None or ratio < worst else worst
    elif passed: extra.append(f"[{a}, {b}] (ratio {ratio.str(3, radius=False)})")
    ok_all &= passed
    print(f"[{a}, {b}]  {u.upper().str(4, radius=False):>12}  {l2.str(4, radius=False):>15}  {lo.str(4, radius=False):>14}  {ratio.str(3, radius=False):>5}  {passed}")
if extra: print("Also passing, not contiguous with the certified range (data points only):", ", ".join(extra))
if end != grid[-1]: print(f"NOT (yet) certified on [{end}, {grid[-1]}]")
if end != grid[0]:
    print(f"CERTIFIED: simple-even for all L in [1.19, {end}] (contiguous Piece III steps, minimum step ratio {worst.str(4, radius=False)});",
          f"with Pieces I-II (code/theorem.py --check): simple-even for all L in (0, {end}]")
else: print("No Piece III step certified")
