# Data

All bounds are Arb balls written as strings `"[midpoint +/- radius]"`. Each window row is
`[lambda0_lower, lambda0_upper, lambda1even_lower, lambda0odd_lower]` for the window `[-L, L]`
(the code certifies the dyadic number `arb(L).upper() >= L`, see paper §4.2).

| File | What it is | Provenance |
|---|---|---|
| `theorem_bounds.json` | **Certificate of record**, 31 windows `0.05 ... 1.19` | Rows `L <= 1.0` copied verbatim from the clean end-to-end rerun (2026-10-24, log `logs/run_rerun.log`). Rows `1.03`, `1.06`, `1.09` copied verbatim from the extension run (2026-10-24/25, log `logs/run_extend.log`; snapshot taken 2026-10-25 after window 1.09 finished). Rows `1.12`, `1.15`, `1.18`, `1.19` copied verbatim from the continued extension run (2026-10-25/26, same log), certified with the upper-triangle panel code. |
| `validate_1.0.json` | Window `L = 1.0` recertified with the sharp-cutoff path (`code/validate10.py`) | Extension run, first step (`logs/run_extend.log`, top). Not part of the gluing. |
| `history/theorem_bounds_original.json` | Bounds from the original window runs (several code revisions, Gauss-Legendre off-by-one) | Superseded; kept for the erratum. |
| `history/theorem_bounds_corrected.json` | The original bounds after the Weyl-shift correction and new upper bounds (`code/gl_correction.py`) | Superseded by the rerun; Table 1 of the paper. |
| `logs/run_rerun.log` | Clean rerun: three stages, the failure at `L = 1.0` with the default margins, the resume with `CERT_MARGIN=1.2 CERT_DDIV=8`, and the gluing | |
| `logs/run_extend.log` | Validation at `L = 1.0`, window `1.03`, the out-of-memory abort at `1.06`, the resume, windows `1.06` and `1.09`, the restart with the upper-triangle code and window `1.12`, then windows `1.15`, `1.18`, `1.19` (cluster assembly by `code/prefetch.py` pipelined with local certification) | Complete run log; the earlier release's log is a prefix of it. |
| `logs/run_gl_correction.log` | Output of `gl_correction.py` (history) | |
| `logs/check.log` | Output of `python code/theorem.py --check` on `theorem_bounds.json` | Ends with `THEOREM CERTIFIED for all L in (0, 1.19]`. |

Notes.

- The extension run kept its own bounds file, which also contains the rows `L <= 1.0`. Those are the same values
  re-serialized, with radii rounded outward by about 1% (for example `3.43e-21` became `3.45e-21`). We release the
  rerun's original strings. Either set passes the gluing.
- The row `L = 1.0` is the rerun's (crude cutoff `T# = 3610`). The released code would recompute `1.0` with the sharp
  path; that result is `validate_1.0.json`.
- Logs were scrubbed for publication: absolute paths are replaced by `<code>/` and `<python>/`, host names by
  `user@host1` ... `user@host3`, and the time-zone abbreviation was removed from timestamps. No number was changed.
