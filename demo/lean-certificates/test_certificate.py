import copy
from fractions import Fraction as Q
import json
from pathlib import Path
import re
import tempfile
import unittest

import checker as c
import check_lean as cl

SOURCE = c.HERE.parent.parent / "data/theorem_bounds.json"


def load(tag):
    return c.parse_data((c.DATA / f"{tag}.lean").read_text())


class CertificateTests(unittest.TestCase):
    """Exact preflight on the committed real certificates, and rejection of tampered ones."""

    def setUp(self):
        self.certs, self.consts = load("W005")

    def reject(self, cert):
        with self.assertRaises(c.Rejected):
            c.check(cert)

    def test_committed_certificates_pass(self):
        for tag in ("W005", "W01"):
            certs, _ = load(tag)
            for name, cert in certs.items():
                with self.subTest(tag=tag, name=name):
                    self.assertTrue(c.check(cert))
                    self.assertEqual(cert["n"], 97)

    def test_emit_parse_roundtrip_is_exact(self):
        for tag in ("W005", "W01"):
            text = (c.DATA / f"{tag}.lean").read_text()
            self.assertEqual(c.emit_data(tag, *c.parse_data(text)), text)

    def test_tampered_matrix_entry(self):
        cert = self.certs["even"]
        cert["C"][40][17] += 2 * cert["K"] + 1
        self.reject(cert)

    def test_wider_interval(self):
        cert = self.certs["odd"]
        cert["rad"] = cert["K"]
        self.reject(cert)

    def test_tampered_factor(self):
        cert = self.certs["deflated"]
        s, vals = cert["Rc"][30]
        cert["Rc"][30] = (s, vals[:-1] + [vals[-1] // 2])
        self.reject(cert)

    def test_mu_too_large(self):
        cert = self.certs["even"]
        cert["mu"] += 1 + cert["M"] - cert["mu"] - cert["n"] * cert["K"]
        self.reject(cert)

    def test_dropping_the_deflation_vector_fails(self):
        cert = self.certs["deflated"]
        cert["v"] = []
        self.reject(cert)

    def test_malformed_rows(self):
        for mutate in (lambda k: k["C"].pop(), lambda k: k["C"][3].append(0), lambda k: k["Rc"].pop(),
                       lambda k: k["Rc"].__setitem__(5, (90, [1] * 10))):
            cert = copy.deepcopy(self.certs["odd"])
            mutate(cert)
            self.reject(cert)

    def test_negative_natural_fields(self):
        cert = self.certs["odd"]
        cert["K"] = -1
        self.reject(cert)

    def test_indefinite_box_cannot_pass(self):
        # A matrix with a negative diagonal entry has no certificate with mu >= 0 for any factor.
        cert = copy.deepcopy(self.certs["even"])
        cert["C"][0][0] = -cert["C"][0][0]
        self.reject(cert)

    def test_pack_unpack(self):
        values = [0, 1, -1, 2 ** 70, -(2 ** 70), 12345]
        wd = c.width(values)
        self.assertEqual(c.unpack(c.pack(values, wd), wd, len(values)), values)
        with self.assertRaises(c.Rejected):
            c.pack([2 ** wd], wd)

    def test_lean_bound_implies_stored_rows(self):
        # min(mu, beta_lo - eD) - eB, the bound Lean derives, is >= the stored upper (hence lower) end.
        rows = c.read_json(SOURCE)
        for tag, L in (("W005", "0.05"), ("W01", "0.1")):
            certs, k = load(tag)
            stored = [c.ball(s) for s in rows[L]]
            for name, col, eD, eB in (("even", 0, "eD_even", "eB_even"), ("deflated", 2, "eD_even", "eB_even"),
                                      ("odd", 3, "eD_odd", "eB_odd")):
                bound = min(Q(certs[name]["mu"], c.SCALE), k["beta_lo"] - k[eD]) - k[eB]
                self.assertGreaterEqual(bound, stored[col][1], (tag, name))

    def test_manifest_matches_data(self):
        manifest = c.read_json(c.HERE / "certificates.json")
        for tag, L in (("W005", "0.05"), ("W01", "0.1")):
            certs, k = load(tag)
            entry = manifest["windows"][L]
            self.assertEqual({n: Q(x) for n, x in entry["mu"].items()},
                             {n: Q(cert["mu"], c.SCALE) for n, cert in certs.items()})
            self.assertEqual({n: Q(x) for n, x in entry["constants"].items()}, k)

    def test_lean_sources_are_clean(self):
        for path in c.HERE.glob("SpectralCertificate/**/*.lean"):
            self.assertIsNone(cl.FORBIDDEN.search(path.read_text()), path.name)


class ParsingTests(unittest.TestCase):
    def test_non_string_numbers_rejected(self):
        for value in (True, 0.0, 1, None, [], {}):
            with self.subTest(value=value), self.assertRaises(c.Rejected):
                c.rational(value)

    def test_invalid_rationals_and_code_rejected(self):
        for value in ("NaN", "Infinity", "1/0", "1/(-2)", "1; sorry", "1e99999", "1" * 513):
            with self.subTest(value=value), self.assertRaises(c.Rejected):
                c.rational(value)

    def test_decimal_is_exact(self):
        self.assertEqual(c.rational("0.1"), Q(1, 10))
        self.assertEqual(c.rational("3.43e-21"), Q(343, 10 ** 23))

    def test_duplicate_json_keys_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.json"
            path.write_text('{"epsilon": "0", "epsilon": "1"}')
            with self.assertRaises(c.Rejected):
                c.read_json(path)

    def test_ball_endpoints_are_outward_exact(self):
        self.assertEqual(c.ball("[0.1 +/- 0.01]"), (Q(9, 100), Q(11, 100)))

    def test_negative_radius(self):
        with self.assertRaises(c.Rejected):
            c.ball("[1 +/- -0.1]")


class StoredRowsTests(unittest.TestCase):
    def test_lean_rows_match_the_table(self):
        # Rows.lean restates lower ends of the stored balls exactly; they must equal data/theorem_bounds.json.
        rows = c.read_json(SOURCE)
        text = (c.HERE / "SpectralCertificate/Rows.lean").read_text()
        S = 10 ** int(re.search(r"def S : ℕ := 10 \^ (\d+)", text)[1])
        for L, tag in (("0.05", "005"), ("0.1", "01")):
            for name, col in (("l1", 2), ("od", 3)):
                with self.subTest(L=L, name=name):
                    value = int(re.search(rf"def {name}_{tag} : ℤ := (\d+)", text)[1])
                    self.assertEqual(Q(value, S), c.ball(rows[L][col])[0])

    def test_window_points_are_the_dyadics_just_above(self):
        points = c.read_json(c.HERE / "certificates.json")["points"]
        self.assertEqual(set(points), {"0.05", "0.1"})
        for L, x in points.items():
            p = c.rational(x)
            self.assertTrue(Q(L) <= p < Q(L) + Q(1, 2 ** 700))
            self.assertEqual(p.denominator & (p.denominator - 1), 0)  # a power of two


class ExporterSmokeTest(unittest.TestCase):
    """The producer must at least import and run its Arb helpers (needs python-flint and the Weil code)."""

    def test_lowest_runs(self):
        try:
            import flint  # noqa: F401
            import export_window as ex
            import sys
            sys.path.insert(0, str(ex.RERUN))
            import certify  # noqa: F401
        except ImportError:
            self.skipTest("python-flint or the Weil code (WEIL_RERUN) is not available")
        C = [[2 * c.SCALE], [0, 3 * c.SCALE]]  # diag(2, 3)
        est, _ = ex.lowest(C, [], 2)
        self.assertAlmostEqual(float(est), 2.0, places=5)  # inverse iteration, 20 steps


class AxiomAuditTests(unittest.TestCase):
    def reports(self, axioms="propext, Classical.choice, Quot.sound"):
        return "\n".join(f"'{name}' depends on axioms: [{axioms}]" for name in cl.EXPECTED)

    def test_standard_foundations(self):
        cl.audit(self.reports())

    def test_rejects_proof_placeholder(self):
        with self.assertRaises(RuntimeError):
            cl.audit(self.reports("propext, sorryAx"))

    def test_rejects_compiler_trust(self):
        with self.assertRaises(RuntimeError):
            cl.audit(self.reports("propext, Lean.ofReduceBool"))

    def test_rejects_custom_axiom(self):
        with self.assertRaises(RuntimeError):
            cl.audit(self.reports("MyAssumedEnclosure"))

    def test_rejects_missing_reports(self):
        with self.assertRaises(RuntimeError):
            cl.audit("")

    def test_no_axioms_report(self):
        cl.audit("\n".join(f"'{name}' does not depend on any axioms" for name in cl.EXPECTED))


if __name__ == "__main__":
    unittest.main()
