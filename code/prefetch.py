"""Assemble head matrices ahead of the certifying run (cache/head_<L>_<p>_<N>.txt), both sectors concurrently.
Run from code/ with the same CLUSTER and CERT_* environment as theorem.py, e.g. `python prefetch.py 1.15 1.18`."""
import sys
from concurrent.futures import ThreadPoolExecutor
import certify

for L in sys.argv[1:]:
    with ThreadPoolExecutor(2) as ex:
        list(ex.map(lambda p: certify.head(L, p, None), (0, 1)))
    print("done", L, flush=True)
