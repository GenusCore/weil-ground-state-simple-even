#!/usr/bin/env python3
"""Exact integer preflight for the Lean certificates, Lean data emission and parsing (stdlib only).

`check` mirrors SpectralCertificate.check in Checker.lean. A Python pass is a preflight, not a proof:
Lean re-evaluates the same conditions in its kernel. `python3 checker.py` re-checks the committed
Lean data files without Arb or the Weil code.
"""
from fractions import Fraction as Q
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
DATA = HERE / "SpectralCertificate" / "Data"
SCALE, HALF = 1 << 64, 1 << 32
KINDS = ("even", "deflated", "odd")
FIELDS = ("n", "C", "Rc", "v", "rad", "M", "K", "mu")


class Rejected(ValueError):
    pass


def rational(value):
    # Exact strings only: no JSON floats, booleans, arbitrary code or Lean terms.
    if not isinstance(value, str) or len(value) > 512:
        raise Rejected("rational must be a string of at most 512 characters")
    if not re.fullmatch(r"[+-]?\d+(?:/\d+|\.\d+(?:[eE][+-]?\d+)?|[eE][+-]?\d+)?", value):
        raise Rejected("invalid rational syntax")
    exponent = re.search(r"[eE]([+-]?\d+)$", value)
    if exponent and abs(int(exponent[1])) > 1000:
        raise Rejected("exponent exceeds limit")
    try:
        return Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise Rejected("invalid rational value") from exc


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise Rejected(f"duplicate JSON field: {key}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=unique)


def ball(value):
    """'[mid +/- rad]' -> exact (mid - rad, mid + rad)."""
    if not isinstance(value, str):
        raise Rejected("stored ball must be a string")
    match = re.fullmatch(r"\[([^\[\]]+) \+/- ([^\[\]]+)\]", value)
    if match is None:
        raise Rejected("unsupported stored ball notation")
    mid, radius = rational(match[1]), rational(match[2])
    if radius < 0:
        raise Rejected("negative ball radius")
    return mid - radius, mid + radius


def entry(C, i, j):
    row = C[max(i, j)] if max(i, j) < len(C) else []
    return row[min(i, j)] if min(i, j) < len(row) else 0


def get(v, i):
    return v[i] if i < len(v) else 0


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def expand(row):
    s, vals = row
    return [0] * s + list(vals)


def resid(c, i, j):
    """Integer residual 2^64 (C + v v^T - M I - R R^T)_ij, as `resid` in Checker.lean."""
    Rc = c["Rc"]
    Ri, Rj = (expand(Rc[i]) if i < len(Rc) else []), (expand(Rc[j]) if j < len(Rc) else [])
    return entry(c["C"], i, j) + get(c["v"], i) * get(c["v"], j) - (c["M"] if i == j else 0) - dot(Ri, Rj)


def residual_bound(c):
    return max(abs(resid(c, i, j)) for i in range(c["n"]) for j in range(i + 1)) + c["rad"]


def check(c):
    """Raise Rejected unless SpectralCertificate.check would return true."""
    if set(c) != set(FIELDS):
        raise Rejected("certificate fields")
    for k in ("n", "rad", "K"):
        if not isinstance(c[k], int) or c[k] < 0:
            raise Rejected(f"{k} must be a natural number")
    if not c["mu"] + c["n"] * c["K"] <= c["M"]:
        raise Rejected("margin: mu + n K > M")
    if len(c["C"]) != c["n"] or len(c["Rc"]) != c["n"]:
        raise Rejected("C and R must have n rows")
    for i, (row, (s, vals)) in enumerate(zip(c["C"], c["Rc"])):
        if len(row) != i + 1 or s < 0 or s + len(vals) > c["n"]:
            raise Rejected(f"malformed row {i}")
    for i in range(c["n"]):
        for j in range(i + 1):
            if abs(resid(c, i, j)) + c["rad"] > c["K"]:
                raise Rejected(f"residual bound fails at ({i}, {j})")
    return True


def compress(R):
    """Dense lower-triangular rows -> (leading zeros, rest up to the diagonal)."""
    out = []
    for i, r in enumerate(R):
        s = next((k for k, x in enumerate(r) if x), i)
        out.append((s, r[s:i + 1]))
    return out


def tag(Ldec):
    return "W" + Ldec.replace(".", "")


# Packed encoding: each list of integers is one hexadecimal natural number, w bits per entry, entry e stored
# as e + 2^(w-1). Lean decodes it with `unpack` (Data.lean) inside the kernel; one big literal elaborates
# in milliseconds, while a list literal with thousands of `Int` numerals takes minutes and gigabytes.
def width(values):
    return max([abs(x) for x in values] + [1]).bit_length() + 2


def pack(values, w):
    off, x = 1 << (w - 1), 0
    for k, e in enumerate(values):
        if not -off <= e < off:
            raise Rejected("value does not fit the packing width")
        x |= (e + off) << (w * k)
    return x


def unpack(x, w, count):
    off, mask = 1 << (w - 1), (1 << w) - 1
    return [((x >> (w * k)) & mask) - off for k in range(count)]


def lean_packed(values):
    w = width(values)
    return f"{1 << w} {1 << (w - 1)}", w, f"0x{pack(values, w):x}"


def emit_data(tagname, certs, consts):
    """One Lean data file per window. The even and deflated certificates share the box."""
    e, d, o = (certs[k] for k in KINDS)
    if d["C"] != e["C"] or d["rad"] != e["rad"] or e["v"] or o["v"] or not (e["n"] == d["n"] == o["n"]):
        raise Rejected("even/deflated must share the box; only the deflated certificate has v")
    n = e["n"]

    def tri(C):
        bo, _, lit = lean_packed([x for r in C for x in r])
        return f"unpackTri {bo} {n} 0\n    {lit}"

    def rows(Rc):
        bo, _, lit = lean_packed([x for _, vals in Rc for x in vals])
        return f"unpackRows {bo} [{', '.join(str(s) for s, _ in Rc)}] 0\n    {lit}"

    def factor(c):
        return f"⟨{rows(c['Rc'])},\n    {c['M']}, {c['K']}, {c['mu']}⟩"

    bo, _, lit = lean_packed(d["v"])
    q = {k: Q(x) for k, x in consts.items()}
    lean_q = {k: f"{x.numerator} / {x.denominator}" for k, x in q.items()}
    return f"""-- Generated by export_window.py from the Arb head-matrix cache. Data only: checked by the kernel.
import SpectralCertificate.Window

namespace SpectralCertificate.{tagname}

def window : Window where
  n := {n}
  Ce := {tri(e["C"])}
  Co := {tri(o["C"])}
  rade := {e["rad"]}
  rado := {o["rad"]}
  v := unpack {bo} {len(d["v"])}\n    {lit}
  fe := {factor(e)}
  fd := {factor(d)}
  fo := {factor(o)}
  betaLo := {lean_q["beta_lo"]}
  eDe := {lean_q["eD_even"]}
  eBe := {lean_q["eB_even"]}
  eDo := {lean_q["eD_odd"]}
  eBo := {lean_q["eB_odd"]}

theorem even_ok : check window.even = true := by decide +kernel

theorem deflated_ok : check window.deflated = true := by decide +kernel

theorem odd_ok : check window.odd = true := by decide +kernel

end SpectralCertificate.{tagname}
"""


def parse_data(text):
    """Decode a generated data file (our own fixed format only) into certificates and constants."""
    def grab(pattern):
        m = re.search(pattern, text, re.M)
        if m is None:
            raise Rejected(f"data file: missing {pattern}")
        return m

    def tri(name):
        m = grab(rf"^  {name} := unpackTri (\d+) (\d+) (\d+) 0\n    0x([0-9a-f]+)$")
        B, n = int(m[1]), int(m[3])
        flat = unpack(int(m[4], 16), B.bit_length() - 1, n * (n + 1) // 2)
        return [flat[i * (i + 1) // 2:(i + 1) * (i + 2) // 2] for i in range(n)]

    def factor(name):
        m = grab(rf"^  {name} := ⟨unpackRows (\d+) (\d+) \[([\d, ]*)\] 0\n    0x([0-9a-f]+),\n    (-?\d+), (\d+), (-?\d+)⟩$")
        B, starts = int(m[1]), [int(s) for s in m[3].split(", ")]
        counts = [i + 1 - s for i, s in enumerate(starts)]
        flat = unpack(int(m[4], 16), B.bit_length() - 1, sum(counts))
        Rc, k = [], 0
        for s, cnt in zip(starts, counts):
            Rc.append((s, flat[k:k + cnt]))
            k += cnt
        return dict(Rc=Rc, M=int(m[5]), K=int(m[6]), mu=int(m[7]))

    n = int(grab(r"^  n := (\d+)$")[1])
    Ce, Co = tri("Ce"), tri("Co")
    m = grab(r"^  v := unpack (\d+) (\d+) (\d+)\n    0x([0-9a-f]+)$")
    v = unpack(int(m[4], 16), int(m[1]).bit_length() - 1, int(m[3]))
    rade, rado = (int(grab(rf"^  {k} := (\d+)$")[1]) for k in ("rade", "rado"))
    certs = {"even": dict(n=n, C=Ce, v=[], rad=rade, **factor("fe")),
             "deflated": dict(n=n, C=Ce, v=v, rad=rade, **factor("fd")),
             "odd": dict(n=n, C=Co, v=[], rad=rado, **factor("fo"))}
    consts = {k: Q(int(grab(rf"^  {lk} := (\d+) / (\d+)$")[1]), int(grab(rf"^  {lk} := (\d+) / (\d+)$")[2]))
              for k, lk in (("beta_lo", "betaLo"), ("eD_even", "eDe"), ("eB_even", "eBe"),
                            ("eD_odd", "eDo"), ("eB_odd", "eBo"))}
    return certs, consts


def main():
    files = sorted(DATA.glob("W*.lean"))
    if not files:
        raise SystemExit("no data files")
    for path in files:
        for name, c in parse_data(path.read_text())[0].items():
            check(c)
            print(f"PASS (python preflight) {path.stem} {name}: n = {c['n']}, mu = {c['mu'] / SCALE:.9f}")
    print("Exact Python preflight passed. This is not Lean verification.")


if __name__ == "__main__":
    try:
        main()
    except Rejected as exc:
        sys.exit(f"REJECTED: {exc}")
