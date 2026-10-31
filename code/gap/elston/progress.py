#!/usr/bin/env python3
"""Job logging with status lines, shared by the native scripts.

Every long-running job writes timestamped lines to logs/<job>.log (and to stdout):

  2026-10-01 14:02:11 | conductor5 sweep | START  args=…
  2026-10-01 14:02:40 | conductor5 sweep | STATUS 3/24 (12.5%)  elapsed 29s  eta 3m23s  | t=tstar x=5 parity=odd
  2026-10-01 14:02:40 | conductor5 sweep | RESULT t=tstar x=5 eps_odd=0.029189
  2026-10-01 14:05:58 | conductor5 sweep | DONE   24/24 in 3m47s

Watch everything:      tail -F logs/*.log
Latest status per job: .venv/bin/python scripts/progress.py   (or: python3 scripts/progress.py)

Usage inside a script:

  from progress import Job
  job = Job("conductor5 sweep", total=len(tasks), args=vars(args))
  for task in tasks:
      ...
      job.step(f"t={t} x={x}")            # advances the counter and writes a STATUS line
      job.result(f"t={t} x={x} eps=…")    # free-form RESULT line
  job.done()
"""

import datetime
import glob
import os
import re
import sys
import time

LOG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs")


def _fmt_dt(seconds):
    seconds = int(max(0, seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h}h{m:02d}m{s:02d}s" if h else (f"{m}m{s:02d}s" if m else f"{s}s")


class Job:
    def __init__(self, name, total=None, args=None, log_path=None, echo=True):
        self.name = name
        self.total = total
        self.count = 0
        self.t0 = time.time()
        self.echo = echo
        os.makedirs(LOG_DIR, exist_ok=True)
        slug = re.sub(r"[^a-zA-Z0-9]+", "-", name).strip("-").lower()
        stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.path = log_path or os.path.join(LOG_DIR, f"{slug}-{stamp}.log")
        self._fh = open(self.path, "a", buffering=1)
        self._write("START ", f"pid={os.getpid()} total={total} log={os.path.relpath(self.path)}" + (f" args={args}" if args else ""))

    def _write(self, kind, text):
        line = f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | {self.name} | {kind} {text}"
        self._fh.write(line + "\n")
        if self.echo:
            print(line, flush=True)

    def log(self, text):
        self._write("INFO  ", text)

    def result(self, text):
        self._write("RESULT", text)

    def warn(self, text):
        self._write("WARN  ", text)

    def step(self, label="", n=1):
        """Advance by n units and write a status line with percent done and ETA."""
        self.count += n
        el = time.time() - self.t0
        if self.total:
            frac = self.count / self.total
            eta = el / frac - el if frac > 0 else 0
            self._write("STATUS", f"{self.count}/{self.total} ({100 * frac:.1f}%)  elapsed {_fmt_dt(el)}  eta {_fmt_dt(eta)}" + (f"  | {label}" if label else ""))
        else:
            self._write("STATUS", f"{self.count} done  elapsed {_fmt_dt(el)}" + (f"  | {label}" if label else ""))

    def sub(self, label, k, n, every=1):
        """Status for an inner loop (does not advance the main counter); writes every `every` iterations."""
        if k % every == 0 or k == n:
            el = time.time() - self.t0
            self._write("STATUS", f"  inner {k}/{n} ({100 * k / n:.1f}%)  elapsed {_fmt_dt(el)}  | {label}")

    def done(self):
        self._write("DONE  ", f"{self.count}/{self.total if self.total else self.count} in {_fmt_dt(time.time() - self.t0)}")
        self._fh.close()


def _last_status():
    rows = []
    for path in sorted(glob.glob(os.path.join(LOG_DIR, "*.log"))):
        with open(path) as fh:
            lines = fh.read().splitlines()
        if not lines:
            continue
        status = next((l for l in reversed(lines) if "| STATUS" in l or "| DONE" in l), lines[-1])
        state = "done" if any("| DONE" in l for l in lines[-3:]) else "running?"
        rows.append((os.path.basename(path), state, status))
    return rows


if __name__ == "__main__":
    rows = _last_status()
    if not rows:
        print(f"no logs in {os.path.relpath(LOG_DIR)}")
        sys.exit(0)
    for name, state, line in rows:
        print(f"{name:<55} {state:<9} {line.split(' | ', 2)[-1]}")
