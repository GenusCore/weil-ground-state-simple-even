"""Cluster worker: sum of certify._partial over panels start, start+step, ... < stop (arg "start:stop:step"),
spread over this machine's cores (process p takes every procs-th panel, so the work stays balanced).
Partial sums are merged by a parallel pairwise tree through files in /dev/shm, so no single process
has to parse every partial. Prints the summed upper triangle, one ball per line."""
import os
import sys
import uuid
from concurrent.futures import ProcessPoolExecutor
from flint import arb
import certify

TMP = f"/dev/shm/weil_{uuid.uuid4().hex}"
DIGITS = int(os.environ.get("CERT_DIGITS", 60))


def work(args):
    flat = certify._partial(args)
    path = f"{TMP}_{uuid.uuid4().hex}"
    with open(path, "w") as f:
        f.write("\n".join(flat))
    return path


def merge(pair):  # streamed line by line, so a merge needs O(1) memory
    a, b = pair
    path = f"{TMP}_{uuid.uuid4().hex}"
    with open(a) as fa, open(b) as fb, open(path, "w") as f:
        for x, y in zip(fa, fb, strict=True):
            f.write((arb(x.rstrip("\n")) + arb(y.rstrip("\n"))).str(DIGITS, radius=True) + "\n")
    os.unlink(a)
    os.unlink(b)
    return path


if __name__ == "__main__":
    L, parity = sys.argv[1], int(sys.argv[2])
    a, b, c = map(int, sys.argv[3].split(":"))
    panels = list(range(a, b, c))
    # half the cores: the coordinator runs the even and odd sectors concurrently on every host.
    # ProcessPoolExecutor (unlike multiprocessing.Pool) raises BrokenProcessPool if a worker dies, e.g. OOM.
    # Memory cap: a worker keeps ~240 N^2 bytes (its O(N^2) arb sums, measured 1.75 GB at N = 2962, plus its
    # /dev/shm file) until the pool closes; half of 80% of MemAvailable per sector. CERT_NODE_PROCS overrides.
    N = certify.Window(L, parity).N
    avail = int(next(l.split()[1] for l in open("/proc/meminfo") if l.startswith("MemAvailable"))) * 1024
    cap = int(os.environ.get("CERT_NODE_PROCS", 0)) or max(1, int(0.4 * avail / (240 * N * N)))
    procs = min(max(1, os.cpu_count() // 2), len(panels), cap)
    print(f"node_partial: N={N} procs={procs}", file=sys.stderr)
    with ProcessPoolExecutor(procs) as pool:
        files = list(pool.map(work, [(L, parity, panels[p::procs]) for p in range(procs)]))
        while len(files) > 1:
            pairs = [(files[i], files[i + 1]) for i in range(0, len(files) - 1, 2)]
            files = list(pool.map(merge, pairs)) + (files[-1:] if len(files) % 2 else [])
    with open(files[0]) as f:
        sys.stdout.write(f.read() + "\n")
    os.unlink(files[0])
