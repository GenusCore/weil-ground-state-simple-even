"""Certified spectral facts about Weil's form for zeta on the window [-L, L] (log scale).

  1. Reproduction (independent code, same method) of Zhu, arXiv:2608.24827 v2, Thm 1.2: Q(f) >= 8.9e-18 ||f||^2, even f, L = 0.8.
  2. "Simple, even ground state" (the spectral hypothesis of Connes-van Suijlekom, arXiv:2511.23257) at several L:
        lambda_0(Q) <= U   (position-space Rayleigh quotient of an explicit trial function)
        lambda_1^even(Q) > U  and  lambda_min^odd(Q) > U   (frequency-side lower bounds).
     Together these give: the lowest eigenvalue is simple and its eigenfunction is even.

Q is the Weil form, as in exploratory/weil_z.py, which checks it against the zeros. For real f with F = f^:
    Q(f) = [pole] + (1/pi) Int_0^inf Psi(t) |F(t)|^2 dt,
    Psi(t) = Re psi(1/4 + it/2) - log pi - sum_{log n < 2L} 2 Lambda(n)/sqrt(n) cos(t log n),
    [pole] = +2F(i/2)^2 for even f and -2|F(i/2)|^2 for odd f (sign flip in the odd sector).
All arithmetic is Arb ball arithmetic. The frequency-cutoff reduction (R), its head/tail block estimate (S) and the
envelope (E) are Zhu's (arXiv:2608.24827, Thm 1.1 and Lemma 3.1), as is the normalization; the proofs below are our
re-derivations, and the code is ours.

(B) Binet: psi(w) = log w - 1/(2w) - Int_0^inf b(x) e^{-wx} dx, b = 1/(e^x-1) - 1/x + 1/2, Re w > 0.
    b(0) = 0 and 0 <= b' <= 1/12. The upper bound reduces to sinh^2 y <= 3y^2/(3-y^2) for y < sqrt3, which holds
    coefficientwise: 2^(2k-1)/(2k)! <= 3^(1-k), with equality at k = 1, 2 and term ratio 4/((2k+1)(2k+2)) <= 1/3.
    Integrating by parts: |psi(w) - log w + 1/(2w)| <= 1/(12 |w| Re w).
(E) Envelope: with w = 1/4 + it/2, Re psi(w) >= log(t/2) - 2/(3t) - 1/(2t^2) >= log(t/2) - 1/t for t >= 3/2.
    The prime comb is >= -A_L, so Psi(t) >= log(T#/2pi) - 1/T# - A_L =: beta for t >= T#.
(R) Plancherel gives (1/pi) Int_0^inf |F|^2 = ||f||^2. Dropping the nonnegative part on [T#, inf):
    Q >= R = [pole] + (1/pi) Int_0^T# (Psi - beta)|F|^2 + beta ||f||^2.  Min-max then gives lambda_k(Q) >= lambda_k(R).
(M) Re psi(s + iy) = psi(s) + sum_n y^2/((n+s)((n+s)^2+y^2)) is increasing in |y|.
(J) |j_n(x)| <= x^n/(2n+1)!!,  i_n(x) <= x^n/(2n+1)!! e^{x^2/2}.
(G) Gauss-Legendre on [-1,1], integrand analytic in the Bernstein ellipse E_rho with |g| <= M:
    error <= (64/15) M rho^(-2n+2)/(rho^2-1) for the n-point rule (Trefethen, ATAP Thm 19.3 states rho^(-2n) for the
    (n+1)-point rule; an earlier revision of this code used rho^(-2n) for n points -- see gl_correction.py).
    Bernstein: a degree-d polynomial bounded by M on [-1,1] is bounded by M rho^d on E_rho.
(S) R = [[A, X], [X^T, D]], D >= (beta - eD) I, ||X|| <= eB  =>  lambda_k(R) >= min(lambda_k(A), beta - eD) - eB.
(D) Deflation: A + v v^T - mu I >= 0  =>  A has at most one eigenvalue below mu (rank-one interlacing).
(P) Position side, for g = f * f~ (even), checked against weil_z.py:
    Q(f) = 4 Int_0^{2L} g cosh(x/2) - g(0)(log pi + gamma + log(1 - e^{-4L}))
           + Int_0^{2L} 2[e^{-2x} g(0) - e^{-x/2} g(x)]/(1 - e^{-2x}) dx - 2 sum Lambda(n)/sqrt(n) g(log n).
    Derivation: psi(s) = -gamma + Int_0^inf (e^{-t} - e^{-st})/(1 - e^{-t}) dt, then Fourier inversion. For x > 2L only
    the g(0) part survives, and it integrates in closed form. For polynomial f, g is a polynomial on [0, 2L].
"""
import os
import subprocess
import sys
import tempfile
import time
from math import e, ceil
from multiprocessing import Pool
from flint import arb, acb, arb_mat, arb_poly, ctx

ctx.prec = int(os.environ.get("CERT_PREC", 768))  # Bessel at x ~ L*T# amplifies input radii by ~e^x: need prec > 1.45 x
NG = 64
B = arb(2) / 5
pi = arb.pi()


def log(x):
    return arb(x).log()


def von_mangoldt(n):
    p = next(d for d in range(2, n + 1) if n % d == 0)
    while n % p == 0:
        n //= p
    return log(p) if n == 1 else None


# (E') For these windows T# is lowered below the crude 2pi e^(A_L + 1/2): the prime comb only nearly reaches A_L at
# much larger t. Psi >= SHARP_BETA on [T#, T2] is certified by interval subdivision (psi_floor); beyond T2 (E) applies.
SHARP_TS = {"1.0": 2450, "1.03": 2450, "1.06": 4000, "1.09": 4000, "1.12": 5550, "1.15": 5550, "1.18": 5550, "1.19": 5550}
SHARP_BETA = arb(3) / 16  # exact dyadic: arb("0.2") parses at low precision and its radius leaks into every entry


def psi_floor(W, lo, hi, h=0.01):
    """Certify Psi(t) >= W.beta for all t in [lo, hi]: Arb enclosures of Psi on cells, bisected where needed.
    No ctx.workprec here: it is global, and head() runs the two parities in threads (a nested save/restore across
    threads left the whole coordinator at the lower precision)."""
    stack = [(arb(lo) + k * arb(h), arb(lo) + (k + 1) * arb(h)) for k in range(ceil((hi - lo) / h))][::-1]
    while stack:
        a, b = stack.pop()
        if W.Psi((a + b) / 2 + arb(0, ((b - a) / 2).upper())) > W.beta:
            continue
        assert (b - a) > arb("1e-9"), f"Psi < beta near t = {a}"
        m = (a + b) / 2
        stack += [(m, b), (a, m)]


class Window:
    def __init__(self, Ldec, parity):
        self.Ldec, self.parity = Ldec, parity
        self.L = L = arb(Ldec).upper()  # exact dyadic >= the decimal
        # frequency cutoff T#: >= 200, and L*T# >= 80 so small windows keep beta well above their eigenvalues
        self.TS = None  # set below, once A_L is known
        self.lam = {}
        n = 2
        while log(n) < 2 * L + 1:
            assert not log(n).overlaps(2 * L)
            if log(n) < 2 * L and von_mangoldt(n):
                self.lam[n] = von_mangoldt(n)
            n += 1
        self.A_L = sum(2 * v / arb(n).sqrt() for n, v in self.lam.items())
        # T#: >= 200, L*T# >= 80 (small windows), and beta = log(T#/2pi) - A_L >= ~0.5 (large windows)
        self.TS = TS = max(200, 5 * ceil(80 / (5 * float(Ldec))),
                           5 * ceil(2 * 3.14159265 * float((self.A_L + arb(0.5)).exp().mid()) / 5))
        self.beta = log(TS) - log(2 * pi) - 1 / arb(TS) - self.A_L
        assert self.beta > 0
        if Ldec in SHARP_TS:  # (E') smaller cutoff; Psi >= SHARP_BETA on [T#, crude T#] is certified by psi_floor in _head
            self.T2, self.TS, self.beta = TS, SHARP_TS[Ldec], SHARP_BETA
            assert self.TS < self.T2 and log(self.T2) - log(2 * pi) - 1 / arb(self.T2) - self.A_L >= SHARP_BETA
        # quadrature nodes per unit panel and extra head modes: the errors N*quad_err and eB must stay << lambda0
        self.NG, extra = (112, 80) if Ldec in SHARP_TS else (NG, 42)
        self.N = ceil(e * float(Ldec) * self.TS / 4) + extra  # head modes; tail orders start past e*L*T#/2

    def order(self, m):
        return 2 * m + self.parity

    def c(self, n):
        return ((2 * n + 1) / (2 * self.L)).sqrt()

    def Psi(self, t):
        return (acb(arb(1) / 4, t / 2).digamma().real - pi.log()
                - sum(2 * v / arb(n).sqrt() * (t * log(n)).cos() for n, v in self.lam.items()))

    def F(self, m, t):  # |F| up to a unit factor, real; Fourier transform of the m-th orthonormal Legendre e_n
        n, x = self.order(m), self.L * t
        j = (pi / (2 * x)).sqrt() * acb(x).bessel_j(arb(2 * n + 1) / 2).real
        return self.c(n) * 2 * self.L * (-1) ** m * j

    def p(self, m):  # pole vector: [pole] = +-(p . a)^2
        n, y = self.order(m), self.L / 2
        return 2 ** arb(.5) * self.c(n) * 2 * self.L * (pi / self.L).sqrt() * y.bessel_i(arb(2 * n + 1) / 2)


def dfact(n):  # (2n+1)!!
    r = arb(1)
    for k in range(1, 2 * n + 2, 2):
        r *= k
    return r


def _partial_lines(args):  # C restricted to the given unit panels [k, k+1]: sum_nodes w (Psi - beta)/pi F F^T
    """Returns the upper triangle (row-major, j >= i) of the panel contribution as ball strings.
    CERT_BESSEL_PREC=b: nodes and Bessel values at b bits via bessel_rec.sph_bessel_hybrid (all orders per node).
    CERT_MATMUL_PREC=m: the N x N product is formed at m bits (rounding goes into the radii)."""
    Ldec, parity, panels = args
    W = Window(Ldec, parity)
    bprec = int(os.environ.get("CERT_BESSEL_PREC", 0))
    if bprec and Ldec in SHARP_TS:  # Bessel at x <= L*T# needs prec > ~1.45 x
        bprec = max(bprec, 512 * ceil(1.75 * float(Ldec) * W.TS / 512))
    mprec = int(os.environ.get("CERT_MATMUL_PREC", ctx.prec))
    if bprec:
        from bessel_rec import sph_bessel_hybrid
        with ctx.workprec(bprec):
            gauss = [arb.legendre_p_root(W.NG, k, weight=True) for k in range(W.NG)]
        nmax = W.order(W.N - 1)
        sign_c = [(-1) ** m * W.c(W.order(m)) * 2 * W.L for m in range(W.N)]
    else:
        gauss = [arb.legendre_p_root(W.NG, k, weight=True) for k in range(W.NG)]

    TINY = arb("1e-100")
    drop = [arb(0)]  # running bound on the error from zeroed entries (added to every entry's radius at the end)

    def column(k, x, w):
        if bprec:
            with ctx.workprec(bprec):
                t = k + (x + 1) / 2
                js = sph_bessel_hybrid(W.L * t, nmax)
            sw, vals = w / 2 * (W.Psi(t) - W.beta) / pi, [sign_c[m] * js[W.order(m)] for m in range(W.N)]
        else:
            t = k + (x + 1) / 2
            sw, vals = w / 2 * (W.Psi(t) - W.beta) / pi, [W.F(m, t) for m in range(W.N)]
        # Zero the negligible tail (|F_m| < 1e-100): wildly mixed magnitudes make FLINT's product ~30x slower.
        # Rigorous: an entry j,k changes by |s|(|F_j||d_k| + |d_j||F~_k|) <= 2 |s| Fmax dmax, summed over nodes.
        cut = W.N
        while cut > 0 and abs(vals[cut - 1]).upper() < TINY:
            cut -= 1
        if cut < W.N:
            dmax = max(abs(v).upper() for v in vals[cut:])
            fmax = max(abs(v).upper() for v in vals)
            drop[0] += 2 * abs(sw).upper() * fmax * dmax
            vals[cut:] = [arb(0)] * (W.N - cut)
        return sw, vals

    # One panel at a time into a running sum of the upper triangle only, kept as row blocks: block b holds rows
    # [bR, bR + R) x columns [bR, N). ~N^2/2 entries per worker instead of three full N x N matrices (C, P, C + P),
    # and half the multiplications. Each entry is the same sum over nodes as the full product G Phi^T.
    R, N = 256, W.N
    C = [None] * ((N + R - 1) // R)
    for k in panels:
        cols = [column(k, x, w) for x, w in gauss]
        with ctx.workprec(mprec):
            for b in range(len(C)):
                r0, r1 = b * R, min(b * R + R, N)
                P = arb_mat([[s * v[m] for s, v in cols] for m in range(r0, r1)]) * \
                    arb_mat([[v[m] for m in range(r0, N)] for _, v in cols])
                C[b] = P if C[b] is None else C[b] + P
    digits = int(os.environ.get("CERT_DIGITS", 240))  # transport precision; rounding is absorbed into the radius
    e = arb(0, drop[0])
    for i in range(N):  # generator: the cluster worker streams these lines to a file
        b, r0 = i // R, i // R * R
        for j in range(i, N):
            yield (C[b][i - r0, j - r0] + e).str(digits, radius=True)


def _partial(args):
    return list(_partial_lines(args))


def _untri(flat, N):  # upper-triangle strings -> full symmetric matrix of arbs
    C = [[None] * N for _ in range(N)]
    it = iter(flat)
    for i in range(N):
        for j in range(i, N):
            C[i][j] = C[j][i] = arb(next(it))
    return C


def assemble_C(Ldec, parity, N, TS, pool):
    """Sum of the panel contributions. With CLUSTER="host1,host2,..." the panels are split across the hosts
    (node_partial.py, all cores of each host); otherwise the local pool does all of them.
    Each host needs a copy of this directory at CERT_REMOTE_DIR (default ~/weil) with a .venv inside it."""
    hosts = [h for h in os.environ.get("CLUSTER", "").split(",") if h]
    if not hosts:
        parts = pool.map(_partial, [(Ldec, parity, range(k, min(k + 5, TS))) for k in range(0, TS, 5)])
        mats = [_untri(P, N) for P in parts]
        return [[sum(M[i][j] for M in mats) for j in range(N)] for i in range(N)]
    # dynamic queue of interleaved panel groups (group g = panels g, g+G, ...): each group mixes cheap and
    # expensive frequencies, and whichever host is free takes the next group
    import hashlib
    import queue
    import threading
    here = os.path.dirname(os.path.abspath(__file__))
    check = "\\n".join(f"{hashlib.sha256(open(os.path.join(here, f), 'rb').read()).hexdigest()}  {f}"
                       for f in ("certify.py", "bessel_rec.py", "node_partial.py"))
    G = 2 * len(hosts)
    todo = queue.Queue()
    for g in range(G):
        todo.put(g)
    C = [[arb(0)] * N for _ in range(N)]
    lock, errors = threading.Lock(), []

    def run(h):
        while True:
            try:
                g = todo.get_nowait()
            except queue.Empty:
                return
            env = " ".join(f"{k}={os.environ[k]}" for k in ("CERT_BESSEL_PREC", "CERT_MATMUL_PREC") if k in os.environ)
            # refuse to run unless the host has byte-identical code (sha256)
            pr = subprocess.run(["ssh", "-o", "BatchMode=yes", h, f"cd {os.environ.get('CERT_REMOTE_DIR', '~/weil')} && printf '{check}\\n' | sha256sum -c --quiet && "
                                 f"CERT_PREC={ctx.prec} CERT_DIGITS=60 {env} "
                                 f".venv/bin/python node_partial.py {Ldec} {parity} {g}:{TS}:{G}"],
                                capture_output=True, text=True)
            vals = pr.stdout.split("\n")
            ntri = N * (N + 1) // 2
            if pr.returncode or len(vals) < ntri:
                errors.append(f"{h} group {g}: {pr.stderr[-300:]}")
                return
            add = _untri(vals[:ntri], N)
            with lock:
                for i in range(N):
                    for j in range(N):
                        C[i][j] += add[i][j]

    threads = [threading.Thread(target=run, args=(h,)) for h in hosts]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors and todo.empty(), errors
    return C


def head(Ldec, parity, pool):
    """Head matrix A of R (with quadrature error folded into the balls) and the tail constants eD, eB.
    Cached in cache/head_<L>_<parity>.txt (ball strings), so a failure after the expensive assembly costs nothing."""
    W = Window(Ldec, parity)
    cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache", f"head_{Ldec}_{parity}_{W.N}.txt")
    if os.path.exists(cache):
        with open(cache) as f:
            vals = f.read().split("\n")
        N = W.N
        eD, eB = arb(vals[0]), arb(vals[1])
        A = [[arb(vals[2 + i * N + j]) for j in range(N)] for i in range(N)]
        return W, A, eD, eB
    W, A, eD, eB = _head(W, Ldec, parity, pool)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    digits = int(ctx.prec * 0.302) + 10
    with open(cache + ".tmp", "w") as f:
        f.write("\n".join([eD.str(digits, radius=True), eB.str(digits, radius=True)]
                          + [x.str(digits, radius=True) for row in A for x in row]))
    os.replace(cache + ".tmp", cache)
    return W, A, eD, eB


def _head(W, Ldec, parity, pool):
    L, N, beta, TS = W.L, W.N, W.beta, W.TS
    if Ldec in SHARP_TS:
        p0 = ctx.prec
        psi_floor(W, TS, W.T2)
        assert ctx.prec == p0
    C = assemble_C(Ldec, parity, N, TS, pool)

    rho = 2 * B + (4 * B * B + 1).sqrt()  # (G) on unit panels, ellipse half-height B
    re_w, abs_w_max = arb(1) / 4 - B / 2, arb(1) / 4 + B / 2 + (arb(TS) + 1) / 2
    psi_e = max(abs_w_max.log(), -re_w.log()) + pi / 2 + 1 / (2 * re_w) + 1 / (12 * re_w * re_w)  # (B)
    comb_e = sum(2 * v / arb(n).sqrt() * (B * log(n)).cosh() for n, v in W.lam.items())
    M = (psi_e + pi.log() + comb_e + beta) * 2 * L * (2 * L * B).exp() / pi
    quad_err = TS * arb(32) / 15 * M * rho ** (-2 * (W.NG - 1)) / (rho ** 2 - 1)

    sign = 1 if parity == 0 else -1
    p = [W.p(m) for m in range(N)]
    A = [[C[i][j] + sign * p[i] * p[j] + (beta if i == j else 0) + arb(0, quad_err) for j in range(N)]
         for i in range(N)]

    X = L * TS  # tails (J), (M), (S)
    n0 = W.order(N)
    sig0 = W.c(n0) * 2 * L * X ** n0 / dfact(n0)
    r = ((2 * n0 + 5) / arb(2 * n0 + 1)).sqrt() * X * X / ((2 * n0 + 3) * (2 * n0 + 5))
    assert r < 1
    sum_sig2 = sig0 ** 2 / (1 - r * r)
    M_psi = (max(abs(acb(arb(1) / 4).digamma().real), abs(acb(arb(1) / 4, arb(TS) / 2).digamma().real))
             + pi.log() + W.A_L + beta)
    y = L / 2
    ptail = 2 ** arb(.5) * W.c(n0) * 2 * L * y ** n0 / dfact(n0) * (y * y / 2).exp() / (1 - y * y).sqrt()
    pH = sum(x * x for x in p).sqrt()
    eD = M_psi / pi * TS * sum_sig2 + (ptail ** 2 if parity else 0)
    eB = M_psi / pi * TS * (2 * L).sqrt() * arb(N).sqrt() * sum_sig2.sqrt() + 2 * pH * ptail
    return W, A, eD, eB


def certified_pd(A, lam):
    if len(A) > 120 and os.path.exists(VCHOL):
        return certified_pd_c(A, lam)
    return certified_pd_py(A, lam)


VCHOL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vchol")


def certified_pd_c(A, lam):
    """Same certificate as certified_pd_py, computed by vchol.c (FLINT, multithreaded)."""
    digits = int(ctx.prec * 0.302) + 10
    with tempfile.NamedTemporaryFile("w", suffix=".mat", delete=False) as f:
        f.write(f"{len(A)}\n{lam.str(digits, radius=True)}\n")
        for row in A:
            f.write("\n".join(x.str(digits, radius=True) for x in row) + "\n")
    prec = os.environ.get("CERT_VCHOL_PREC", str(ctx.prec))
    out = subprocess.run([VCHOL, f.name, prec, str(os.cpu_count()), os.environ.get("CERT_DDIV", "100")], capture_output=True, text=True).stdout
    os.unlink(f.name)
    return out.startswith("PD")


def certified_pd_py(A, lam):
    """A - lam*I >= 0. Approximate Cholesky of A - (lam+d)I on midpoints (no interval blowup), then a rigorous
    residual E = A - (lam+d)I - R R^T: ||E||_F < d  =>  A - lam*I = R R^T + E + dI >= 0."""
    n, d = len(A), lam / 100
    R = [[arb(0)] * n for _ in range(n)]
    for j in range(n):
        s = (A[j][j] - lam - d - sum(R[j][k] ** 2 for k in range(j))).mid()
        if not s > 0:
            return False
        R[j][j] = s.sqrt().mid()
        for i in range(j + 1, n):
            R[i][j] = ((A[i][j] - sum(R[i][k] * R[j][k] for k in range(j))) / R[j][j]).mid()
    Rm = arb_mat(R)
    E = arb_mat(A) - arb_mat([[lam + d if i == j else 0 for j in range(n)] for i in range(n)]) - Rm * Rm.transpose()
    return sum(abs(E[i, j]).upper() ** 2 for i in range(n) for j in range(n)).sqrt() < d


def lowest(A, iters=20):
    """Approximate lowest eigenpair by inverse iteration on midpoints (only used to choose what to certify)."""
    n = len(A)
    veig = os.path.join(os.path.dirname(os.path.abspath(__file__)), "veig")
    if n > 300 and os.path.exists(veig):  # one LU in C instead of `iters` full solves
        with tempfile.NamedTemporaryFile("w", suffix=".mat", delete=False) as f:
            f.write(f"{n}\n" + "\n".join(x.mid().str(110, radius=False) for row in A for x in row) + "\n")
        out = subprocess.run([veig, f.name, "384", str(iters)], capture_output=True, text=True).stdout.split()
        os.unlink(f.name)
        v = [arb(x).mid() for x in out[1:n + 1]]
        return sum(v[i] * A[i][j] * v[j] for i in range(n) for j in range(n)).mid(), v
    with ctx.workprec(320):
        Am = arb_mat([[x.mid() for x in row] for row in A])
        x = arb_mat([[1] for _ in range(n)])
        for _ in range(iters):
            x = Am.solve(x)
            nrm = sum(x[i, 0] ** 2 for i in range(n)).sqrt()
            x = arb_mat([[(x[i, 0] / nrm).mid()] for i in range(n)])
    v = [x[i, 0] for i in range(n)]
    return sum(v[i] * A[i][j] * v[j] for i in range(n) for j in range(n)).mid(), v


def certify_lower(A, est, eD, eB, beta):
    mu = est / arb(os.environ.get("CERT_MARGIN", "1.05"))  # certify a hair below the estimate
    assert certified_pd(A, mu), "could not certify"
    return min(mu, beta - eD) - eB


def rayleigh_upper(W, v, n_out=450, rho=arb(3), K=2048):
    """(P) + (G): rigorous upper bound on Q(f)/||f||^2 for f = sum v_m e_m (even, polynomial)."""
    with ctx.workprec(1536):
        L = W.L
        x = arb_poly([0, 1])
        P = [arb_poly([1]), x]
        for k in range(1, 2 * len(v)):
            P.append(((2 * k + 1) * x * P[k] - k * P[k - 1]) * (arb(1) / (k + 1)))
        f = sum((v[m] * W.c(2 * m) * P[2 * m] for m in range(len(v))), arb_poly([0]))  # f(Ls) as poly in s
        D = 2 * len(v) - 2
        inner = [arb.legendre_p_root(D + 1, k, weight=True) for k in range(D + 1)]

        def g(s):  # g(L s) = L Int_{-1}^{1-s} f(y) f(y+s) dy, exact (polynomial of degree 2D)
            h = (2 - s) / 2
            return L * h * sum(w * f(-1 + h * (y + 1)) * f(-1 + h * (y + 1) + s) for y, w in inner)

        norm2 = sum(c * c for c in v)
        g0 = g(arb(0))

        def H(s, gs):  # (P) integrand on [0, 2] in the variable s = x/L
            q = (-2 * L * s).exp()
            return 4 * L * (L * s / 2).cosh() * gs + 2 * L * (q * g0 - (-L * s / 2).exp() * gs) / (1 - q)

        outer = [arb.legendre_p_root(n_out, k, weight=True) for k in range(n_out)]
        integral = sum(w * H(1 + y, g(1 + y)) for y, w in outer)

        # (G) for H on the ellipse with foci 0, 2: |g| <= rho^(2D+1) ||f||^2 (Bernstein; |g| <= ||f||^2 on [0,2])
        a_, b_ = (rho + 1 / rho) / 2, (rho - 1 / rho) / 2
        Gb = rho ** (2 * D + 1) * norm2
        Mh = arb(0)
        for k in range(K):
            th = arb(2 * k + 1) * pi / K + arb(0, pi / K)  # covers the k-th arc
            z = acb(1 + a_ * th.cos(), b_ * th.sin())
            den = abs(1 - (-2 * L * z).exp()).lower()
            assert den > 0
            num = abs((-2 * L * z).exp()) * abs(g0) + abs((-L * z / 2).exp()) * Gb
            val = (4 * L * abs((L * z / 2).cosh()) * Gb + 2 * L * num / den).upper()
            Mh = val if val > Mh else Mh
        quad = arb(64) / 15 * Mh * rho ** (-2 * (n_out - 1)) / (rho ** 2 - 1)

        gamma = arb.const_euler()
        Q = (integral + arb(0, quad) - g0 * (pi.log() + gamma + (1 - (-4 * L).exp()).log())
             - 2 * sum(lv / arb(n).sqrt() * g(log(n) / L) for n, lv in W.lam.items()))
        return (Q / norm2).upper()


def window_bounds(Ldec, pool):
    """Certified (lam0_lo, lam0_up, lam1even_lo, lam0odd_lo) for the window [-L, L]."""
    if os.environ.get("CLUSTER"):  # keep the cluster queue full: assemble both sectors concurrently
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(2) as ex:
            fe, fo = ex.submit(head, Ldec, 0, pool), ex.submit(head, Ldec, 1, pool)
            (We, Ae, eDe, eBe), (Wo, Ao, eDo, eBo) = fe.result(), fo.result()
    else:
        We, Ae, eDe, eBe = head(Ldec, 0, pool)
        Wo, Ao, eDo, eBo = head(Ldec, 1, pool)
    l0_est, v = lowest(Ae)
    lam0_lo = certify_lower(Ae, l0_est, eDe, eBe, We.beta)
    n = len(v)
    Ad = [[Ae[i][j] + v[i] * v[j] for j in range(n)] for i in range(n)]  # (D)
    lam1_lo = certify_lower(Ad, lowest(Ad)[0], eDe, eBe, We.beta)
    odd_lo = certify_lower(Ao, lowest(Ao)[0], eDo, eBo, Wo.beta)
    U = rayleigh_upper(We, v[:150])  # the ground state lives on low orders; a truncated trial is still a valid trial
    assert lam0_lo <= U
    if Ldec == "0.8":
        assert certified_pd(Ae, arb("8.9e-18")) and min(arb("8.9e-18"), We.beta - eDe) - eBe > 0  # Zhu Thm 1.2
    return lam0_lo, U, lam1_lo, odd_lo


if __name__ == "__main__":
    t0 = time.time()
    with Pool(20) as pool:
        print(" L     lambda_0 >=   lambda_0 <=   lambda_1^even >=  lambda_min^odd >=   simple & even?")
        for Ldec in sys.argv[1:] or ["0.4", "0.5", "0.6", "0.7", "0.8"]:
            lo, U, l1, od = window_bounds(Ldec, pool)
            print(f"{Ldec:4}  {lo.str(4):>12}  {U.str(4):>12}  {l1.str(4):>14}  {od.str(4):>16}   {l1 > U and od > U}"
                  f"   ({time.time() - t0:.0f}s)", flush=True)
