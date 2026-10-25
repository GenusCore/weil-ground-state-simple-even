"""Piece I: Lipschitz control of the rescaled Weil form for small windows (no primes: L <= log(2)/2).

Rescale f(u) = L^{-1/2} f1(u/L), ||f1|| = 1 on [-1, 1]. Then Q(f) = log(1/L) + R_L(f1), with
    R_L(f1) = P_L(f1) + (1/pi) Int_0^inf m_L(tau) |F1(tau)|^2 dtau,
    m_L(tau) = Re psi(1/4 + i tau/(2L)) - log pi + log L,
    P_L = +2L (Int f1 e^{-Ls/2})^2 (even f1),   -2L (Int f1 sinh(Ls/2))^2 (odd f1).
d m_L / dL = (1/L) h(tau/(2L)),  h(u) = 1 - Re(i u psi'(1/4 + iu)).
Since |F1|^2 <= (Int|f1|)^2 <= 2, the multiplier part of dR/dL is bounded by (4/pi) Int_0^inf |h|.

Bound on h. psi'(w) = sum_n (n+w)^-2 gives h(u) = 1 - sum_{n>=0} phi(n + 1/4), with phi(a) = 2a u^2/(a^2+u^2)^2
and Int_0^inf phi = 1. The midpoint rule on cells of width 1 (Peano kernel <= 1/8), together with
Int_R |phi''| = 6/u^2 and |Int_{-1/4}^0 phi| <= 1/(16u^2), gives |h(u)| <= (3/4 + 1/16)/u^2 = 0.8125/u^2.

The eigenvalues of R_L in each parity sector are therefore Lipschitz in L with constant C_sector (min-max over the common
core C_c^inf). Gaps between eigenvalues of Q_L equal gaps of R_L.
"""
from flint import arb, acb

pi = arb.pi()


def h_sup(a, b):  # sup |h| on [a, b] by ball evaluation
    u = arb((a + b) / 2, (b - a) / 2)
    return abs(1 - (acb(0, u) * acb(arb(1) / 4, u).polygamma(1)).real).upper()


def int_abs_h(U0=40, per_unit=200):
    n, step = U0 * per_unit, arb(1) / per_unit
    head = sum(step * h_sup(k * step, (k + 1) * step) for k in range(n))
    return head + arb("0.8125") / U0


def constants(L0):
    L0 = arb(L0)
    I = int_abs_h()
    mult = 4 / pi * I
    J2, Jp = 2 * L0.sinh() / L0, (2 * L0.exp() / 3).sqrt() / 2  # |J|^2, |J'| bounds for the even pole
    C_e = mult + 2 * J2 + 4 * L0 * J2.sqrt() * Jp
    K2, Kp = 2 * (L0 / 2).sinh() ** 2, (L0 / 2).cosh() * (arb(2) / 3).sqrt() / 2  # odd pole
    C_o = mult + 2 * K2 + 4 * L0 * K2.sqrt() * Kp
    return I, C_e.upper(), C_o.upper()


if __name__ == "__main__":
    I, C_e, C_o = constants("0.05")
    print(f"Int|h| <= {I.upper().str(6)}   C_even <= {C_e.str(6)}   C_odd <= {C_o.str(6)}  (valid for L <= 0.05)")
    assert abs(h_sup(0, arb(1) / 1000) - 1) < arb("0.01")  # h(0) = 1
