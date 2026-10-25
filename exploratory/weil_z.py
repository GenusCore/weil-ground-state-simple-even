"""Weil's quadratic form for zeta: the Z analogue of the Hodge-index check in weil.py.

For real f supported in [-a, a] (log scale, x = log of a multiplicative variable),
    Q(f) = W(f * f~),   W = Weil explicit-formula functional.
Explicit formula (prime side):
    W(g) = gh(i/2) + gh(-i/2) - g(0) log(pi) + (1/2pi) Int gh(r) Re psi(1/4 + ir/2) dr
           - sum_n Lambda(n)/sqrt(n) [g(log n) + g(-log n)],     gh(t) = Int g(u) e^{itu} du
Zero side: W(g) = sum over zeros rho = 1/2 + i*gamma of gh(gamma).
Weil: RH <=> Q(f) >= 0 for every f, every a. Prime p^k enters once 2a > k log p.

Q is computed from the PRIME side only. The zero side (first zeros from mpmath) is used
purely as a test of the implementation. Finite basis => Rayleigh-Ritz: our lambda_min is an
UPPER bound on the true minimum. It can refute positivity, never certify it.
"""
import numpy as np
from numpy.polynomial import legendre as L
from pathlib import Path

DEG, NODES, NZEROS = 14, 800, 400
T, TW = L.leggauss(NODES)


def digamma(z):  # complex, Re z > 0: shift up, then asymptotic series
    z, acc = np.asarray(z, complex), 0
    for _ in range(10):
        acc, z = acc - 1 / z, z + 1
    w = 1 / z ** 2
    return acc + np.log(z) - .5 / z - w * (1 / 12 - w * (1 / 120 - w * (1 / 252 - w * (1 / 240 - w / 132))))


def von_mangoldt(N):
    lam = np.zeros(N + 1)
    for p in range(2, N + 1):
        if all(p % d for d in range(2, int(p ** .5) + 1)):
            pk = p
            while pk <= N:
                lam[pk], pk = np.log(p), pk * p
    return lam


class Basis:
    """phi_n(x) = (1 - (x/a)^2)^3 P_n(x/a) on [-a, a]; C^2, so its Fourier transform decays like r^-4."""

    def __init__(self, a, deg=DEG):
        self.a, self.deg = a, deg
        self.x, self.w = a * T, a * TW
        self.Phi = self(self.x)
        self.M = self.Phi * self.w @ self.Phi.T  # <phi_j, phi_k>

    def __call__(self, x):
        t = np.asarray(x, float) / self.a
        inside = np.abs(t) <= 1
        return np.array([np.where(inside, (1 - t * t) ** 3 * L.legval(t, np.eye(self.deg)[n]), 0)
                         for n in range(self.deg)])

    def ft(self, t):  # phi_hat(t) = Int phi(x) e^{itx} dx, t may be complex
        return (self.Phi * self.w) @ np.exp(1j * np.outer(self.x, np.atleast_1d(t)))

    def corr(self, u):  # g_jk(u) = Int phi_j(x) phi_k(x - u) dx
        return self.Phi * self.w @ self(self.x - u).T


def ghat(B, t):  # ghat_jk(t) = phi_hat_j(t) phi_hat_k(-t), summed over the points t
    t = np.atleast_1d(t)
    return (B.ft(t) @ B.ft(-t).T).real


def archimedean(B, rmax=None, panel=.5, k=16):
    # ponytail: quadrature ft is only exact while r*a << NODES, so cut there; tail ~ (ar)^-8
    rmax = rmax or 0.5 * NODES / B.a
    gx, gw = L.leggauss(k)
    A = np.zeros((B.deg, B.deg))
    for lo in np.arange(0, rmax, 400.):
        edges = np.arange(lo, lo + 400., panel)
        r = (edges[:, None] + panel / 2 * (gx + 1)).ravel()
        wr = np.tile(panel / 2 * gw, len(edges))
        F = B.ft(r)
        A += ((F * (wr * digamma(.25 + .5j * r).real)) @ F.conj().T).real
    return A / np.pi  # integrand is even in r: (1/2pi) Int_R = (1/pi) Int_0^inf


def Q_primes(B, weight=1.):
    Q = ghat(B, [.5j, -.5j]) - np.log(np.pi) * B.M + archimedean(B)
    lam = von_mangoldt(int(np.exp(2 * B.a)) + 1)
    for n in np.nonzero(lam)[0]:
        if np.log(n) < 2 * B.a:
            G = B.corr(np.log(n))
            Q -= weight * lam[n] / np.sqrt(n) * (G + G.T)
    return Q


def Q_zeros(B, gammas):
    F = B.ft(gammas)
    return 2 * (F @ F.conj().T).real


def lam_min(B, Q):  # min of Q(f)/||f||^2 over the basis span
    Linv = np.linalg.inv(np.linalg.cholesky(B.M))
    return np.linalg.eigvalsh(Linv @ Q @ Linv.T)[0]


def zeros(n=NZEROS, cache=Path(__file__).with_name("zeros.txt")):
    if not cache.exists():
        import mpmath
        cache.write_text("\n".join(str(float(mpmath.zetazero(k).imag)) for k in range(1, n + 1)))
    return np.array(cache.read_text().split(), float)[:n]


if __name__ == "__main__":
    assert abs(digamma(1).real + 0.5772156649015329) < 1e-13

    # Test 1: prime side == zero side (validates every term of the explicit formula)
    gam = zeros()
    for a in (0.2, 0.5, 0.9):  # no primes / primes 2,(4) / up to 5
        B = Basis(a)
        Qp, Qz = Q_primes(B), Q_zeros(B, gam)
        err = np.abs(Qp - Qz).max()
        print(f"test explicit formula a={a}: max |prime side - zero side| = {err:.1e}")
        assert err < 1e-10, err

    # Test 2: tail of the archimedean integral is negligible
    for a in (0.2, 1.2):
        B = Basis(a)
        assert np.abs(archimedean(B, 0.25 * NODES / a) - archimedean(B)).max() < 1e-12

    # Sweep: smallest eigenvalue as support grows; mark prime powers entering
    lam = von_mangoldt(20)
    print("\n   a   prime powers in play    lam_min(prime side)  lam_min(zero side)  lam_min(control)")
    ctrl = [20 - .25j, 20 + .25j, -20 - .25j, -20 + .25j]  # extra off-line zero quadruple, beta = 1/2 +- 1/4
    for a in np.arange(.1, 1.21, .1):
        B = Basis(a)
        Qp = Q_primes(B)
        pp = [int(n) for n in np.nonzero(lam)[0] if np.log(n) < 2 * a]
        print(f"{a:5.2f}  {str(pp):24} {lam_min(B, Qp):18.3e} {lam_min(B, Q_zeros(B, gam)):20.3e}"
              f" {lam_min(B, Qp + ghat(B, ctrl)):18.3e}")
