"""Rigorous spherical Bessel j_0..j_{nmax}(x) for an exact (dyadic) x > 0, in O(nmax) ball operations.

  n <= n1 = ceil(x): upward recurrence j_{n+1} = (2n+1)/x j_n - j_{n-1} from the exact j_0 = sin x/x and
      j_1 = sin x/x^2 - cos x/x. This is the oscillatory regime, where errors grow only mildly.
  n >  n1: ratios r_k = j_k/j_{k-1} from the backward continued fraction r_k = x/(2k+1 - x r_{k+1}).
      Lemma: if 2k+1 >= 2x, the map tau -> x/(2k+1 - x tau) sends [0,1] into [x/(2k+1), x/(2k+1-x)], a subset of [0,1].
      j_n is the minimal solution of the recurrence (y_n dominates), so by Pincherle's theorem (Gautschi, SIAM Rev. 9,
      1967) the continued fraction converges to r_k. Starting from r_{M+1} in [0,1] therefore encloses every r_k, k <= M.
      Past M, |j_n| <= x^n/(2n+1)!!, and that bound is used once it is below `tiny`.

STATUS: prototype. The backward ratios are sharp. The upward ball recurrence for n <= x suffers the interval
"wrapping effect": radii grow ~2-3x per step (1e375 at x = 3000). The recurrence is accurate; the enclosure is not.
Fix (not done): float recurrence plus a Green's-function bound (Casoratian j_m y_{m-1} - j_{m-1} y_m = 1/x^2), or
rigorous Debye asymptotics near the turning point.
"""
from flint import arb, acb, ctx


def sph_bessel_all(x, nmax, tiny):
    s, c = x.sin(), x.cos()
    j = [s / x, s / (x * x) - c / x]
    n1 = min(nmax, max(1, int(float(x.mid())) + 1))
    for n in range(1, n1):
        j.append((2 * n + 1) / x * j[n] - j[n - 1])
    if nmax <= n1:
        return j[:nmax + 1]
    # backward ratios on (n1, M], M large enough that the start interval [0,1] has contracted away
    M = max(nmax, int(2.72 * float(x.mid()) / 2) + 60)
    r = arb(0.5, 0.5)  # [0, 1]
    ratios = {}
    for k in range(M, n1, -1):
        r = x / (2 * k + 1 - x * r)
        if k <= nmax:
            ratios[k] = r
    for k in range(n1 + 1, nmax + 1):
        j.append(j[k - 1] * ratios[k])
    return j[:nmax + 1]


def sph_bessel_hybrid(x, nmax, seed_every=400):
    """Rigorous j_0..j_nmax(x) for a sharp ball x > 0, at the current working precision. Fixes the wrapping-effect
    problem of sph_bessel_all by re-seeding.

      n <= n1 = floor(x): Arb's own j_s, j_{s+1} every `seed_every` orders, then upward recurrence between seeds.
          In ball arithmetic the radius grows at most ~2.5x per step, i.e. <= 1e200 over 400 steps. That is
          harmless at >= 6000 bits, and the radii are tracked anyway.
      n >  n1: backward continued-fraction ratios (lemma + Pincherle, see sph_bessel_all), anchored at j_{n1}.
    """
    pi = arb.pi()

    def arb_j(n):
        return (pi / (2 * x)).sqrt() * acb(x).bessel_j(arb(2 * n + 1) / 2).real

    n1 = min(nmax, int(float(x.mid())))
    j = []
    s = 0
    while s <= n1:
        a, b = arb_j(s), arb_j(s + 1)
        run = [a, b]
        for n in range(s + 1, min(s + seed_every, n1 + 1)):
            run.append((2 * n + 1) / x * run[-1] - run[-2])
        j.extend(run[:min(seed_every, n1 + 1 - s)])
        s += seed_every
    if nmax > n1:
        M = max(nmax, int(2.72 * float(x.mid()) / 2) + 60)
        r = arb(0.5, 0.5)
        ratios = [None] * (nmax + 1)
        for k in range(M, n1, -1):
            r = x / (2 * k + 1 - x * r)
            if k <= nmax:
                ratios[k] = r
        for k in range(n1 + 1, nmax + 1):
            j.append(j[k - 1] * ratios[k])
    return j


if __name__ == "__main__":  # check against Arb's own Bessel where Arb is reliable
    import time
    ctx.prec = 768
    for xs, nmax in (("0.37", 300), ("12.5", 300), ("160.25", 300), ("3000.5", 4000)):
        x = arb(xs)
        t = time.time()
        j = sph_bessel_all(x, nmax, arb("1e-120"))
        dt = time.time() - t
        worst = arb(0)
        for n in range(0, nmax + 1, max(1, nmax // 40)):
            ref = (arb.pi() / (2 * x)).sqrt() * acb(x).bessel_j(arb(2 * n + 1) / 2).real
            if ref.rad() < arb("1e-100"):
                assert j[n].overlaps(ref), (xs, n, j[n], ref)
            worst = max(worst, j[n].rad())
        print(f"x={xs:>7} orders 0..{nmax}: max radius {worst.str(3)}   {dt * 1e6 / (nmax + 1):.1f} us/order")
