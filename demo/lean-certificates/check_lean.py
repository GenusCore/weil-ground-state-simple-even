#!/usr/bin/env python3
"""Build the Lean library, audit the axioms of the theorems, and ask Lean to reject tampered data.

Requires Lean/Lake (lean-toolchain) and the Mathlib dependency (`lake exe cache get`). Exits nonzero on any
failure, including a missing toolchain. Never turns a Python pass into a Lean pass.
"""
import copy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

import checker as c

EXPECTED = [
    "SpectralCertificate.check_sound", "SpectralCertificate.check_sound_perp",
    "SpectralCertificate.window_sound",
    *(f"SpectralCertificate.{t}.{k}_ok" for t in ("W005", "W01") for k in c.KINDS),
    "SpectralCertificate.W005_rows", "SpectralCertificate.W01_rows",
]
ALLOWED = {"propext", "Classical.choice", "Quot.sound"}
FORBIDDEN = re.compile(r"\b(sorry|admit|native_decide|implemented_by|extern|unsafe|axiom|ofReduceBool)\b")


def run(args, check=True):
    t0 = time.time()
    result = subprocess.run(args, cwd=c.HERE, text=True, capture_output=True)
    out = result.stdout + result.stderr
    if check and result.returncode:
        print(out)
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(args)}")
    return out, result.returncode, time.time() - t0


def audit(text):
    found = set()
    for match in re.finditer(r"'([^']+)' depends on axioms:\s*\[([^\]]*)\]", text):
        name, raw = match.groups()
        axioms = set(filter(None, (x.strip() for x in raw.split(","))))
        if name in EXPECTED:
            if not axioms <= ALLOWED:
                raise RuntimeError(f"unexpected axioms for {name}: {sorted(axioms - ALLOWED)}")
            found.add(name)
    for name in EXPECTED:
        if f"'{name}' does not depend on any axioms" in text:
            found.add(name)
    if found != set(EXPECTED):
        raise RuntimeError(f"missing axiom reports: {sorted(set(EXPECTED) - found)}")


def lean_file(directory, name, text):
    path = Path(directory) / f"{name}.lean"
    path.write_text(text)
    return path


def only(text, theorem):
    """Keep one of the three kernel checks in a data file (faster negative tests)."""
    return re.sub(rf"^theorem (?!{theorem}\b)\w+ : check .*$\n", "", text, flags=re.M)


def negative_examples(directory):
    data = (c.DATA / "W005.lean").read_text()
    certs, consts = c.parse_data(data)
    cases = []

    def tampered(label, kind, mutate):
        bad = copy.deepcopy(certs)
        mutate(bad)
        if kind == "deflated":  # keep the shared box consistent so only the intended change is tested
            bad["even"]["C"], bad["even"]["rad"] = bad["deflated"]["C"], bad["deflated"]["rad"]
        cases.append((label, only(c.emit_data("Tampered", bad, consts), f"{kind}_ok"),
                      f"proved that the proposition\n  check window.{kind} = true\nis false"))

    def bump_entry(b):
        b["even"]["C"][40][17] += 2 * b["even"]["K"] + 1  # one center entry moved past the residual budget
    tampered("matrix-entry", "even", bump_entry)

    def widen(b):  # K = max |residual| + rad exactly, so one more unit of radius exceeds the budget
        assert b["odd"]["K"] == c.residual_bound(b["odd"])
        b["odd"]["rad"] += 1
    tampered("interval-radius+1", "odd", widen)

    def factor(b):
        s, vals = b["deflated"]["Rc"][30]
        b["deflated"]["Rc"][30] = (s, vals[:-1] + [vals[-1] // 2])  # one factor entry (a diagonal one) halved
    tampered("factor-entry", "deflated", factor)

    def raise_mu(b):  # one unit (2^-64) above the largest mu the margin allows
        b["even"]["mu"] = b["even"]["M"] - b["even"]["n"] * b["even"]["K"] + 1
    tampered("mu+1", "even", raise_mu)

    # Boundary acceptance: the largest admissible mu (exactly M - n K) must still pass, so the two
    # rejections above are one unit past the edge, not far outside it.
    edge = copy.deepcopy(certs)
    edge["even"]["mu"] = edge["even"]["M"] - edge["even"]["n"] * edge["even"]["K"]
    path = lean_file(directory, "Boundary_mu", only(c.emit_data("Boundary", edge, consts), "even_ok"))
    out, code, secs = run(["lake", "env", "lean", str(path)], check=False)
    if code:
        raise RuntimeError(f"Lean rejected the boundary mu = M - n K:\n{out[-2000:]}")
    print(f"PASS: Lean accepts the boundary mu = M - n K ({secs:.1f} s)")
    results = [{"case": "accept mu = M - n K", "seconds": round(secs, 1)}]
    for label, text, expected in cases:
        path = lean_file(directory, f"Negative_{label.replace('-', '_')}", text)
        out, code, secs = run(["lake", "env", "lean", str(path)], check=False)
        if code == 0 or expected not in out:
            raise RuntimeError(f"Lean failed to reject {label}:\n{out[-2000:]}")
        print(f"PASS: Lean rejected tampered {label} ({secs:.1f} s)")
        results.append({"case": label, "seconds": round(secs, 1)})
    return results


def main():
    output = c.HERE / "validation-local.json"
    output.unlink(missing_ok=True)  # never leave a stale success report behind
    if shutil.which("lake") is None:
        raise SystemExit("Lean verification NOT RUN: lake is not installed. See README.md.")
    sources = sorted(c.HERE.glob("SpectralCertificate/**/*.lean")) + [c.HERE / "SpectralCertificate.lean"]
    for path in sources:
        if FORBIDDEN.search(path.read_text()):
            raise RuntimeError(f"forbidden keyword in {path.name}")
    _, _, build_secs = run(["lake", "build"])
    print(f"PASS: lake build ({build_secs:.1f} s)")
    with tempfile.TemporaryDirectory(prefix="lean-check-", dir=c.HERE) as directory:
        probe = lean_file(directory, "Audit", "import SpectralCertificate\n" +
                          "".join(f"#print axioms {name}\n" for name in EXPECTED))
        out, _, _ = run(["lake", "env", "lean", str(probe)])
        audit(out)
        print(f"PASS: axiom audit of {len(EXPECTED)} theorems (only {', '.join(sorted(ALLOWED))})")
        negatives = negative_examples(directory)
    hashes = {str(p.relative_to(c.HERE)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    manifest = c.HERE / "lake-manifest.json"
    evidence = {"status": "lean-checked", "utc": datetime.now(timezone.utc).isoformat(),
                "toolchain": (c.HERE / "lean-toolchain").read_text().strip(),
                "lean_version": run(["lake", "env", "lean", "--version"])[0].strip(),
                "lake_build_seconds": round(build_secs, 1),
                "axiom_audit": EXPECTED, "allowed_axioms": sorted(ALLOWED),
                "negative_tests": negatives, "lean_source_sha256": hashes,
                "lake_manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest()}
    output.write_text(json.dumps(evidence, indent=2) + "\n")
    print("PASS: Lean build, axiom audit and negative tests")


if __name__ == "__main__":
    main()
