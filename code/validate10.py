"""Validation of the sharp-cutoff code path (certify.SHARP_TS, NG=112): recertify L = 1.0 with T# = 2450 (N = 1745)
and compare with the clean rerun (crude T# = 3610, N = 2496). Writes validate_1.0.json; touches nothing else."""
import json
from multiprocessing import Pool
from flint import arb
import certify
with Pool(20) as pool:
    r = certify.window_bounds("1.0", pool)
json.dump([x.str(20, radius=True) for x in r], open("validate_1.0.json", "w"))
old = [arb(x) for x in json.load(open("theorem_bounds.json"))["1.0"]]
for name, a, b in zip(["lam0_lo", "lam0_up", "lam1even_lo", "lam0odd_lo"], r, old):
    print(f"{name:12} sharp {a.str(6)}   rerun {b.str(6)}   ratio {(a / b).str(4)}")
lo, U, l1, od = r
assert lo <= U and l1 > U and od > U, "validation FAILED"
assert lo <= old[1] and old[0] <= U, "brackets inconsistent with the rerun"
print("VALIDATION OK")
