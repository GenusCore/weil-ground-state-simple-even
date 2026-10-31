"""Convergence sweep with the fast recurrence: python fr_sweep.py out.jsonl "K1,K2" "L1,L2" prec. One process per (K, L)."""
import sys, json, time
from multiprocessing import Pool
def job(a):
    K, L, prec = a
    import hp, fast_recurrence as fr
    from flint import ctx, arb
    import os; ctx.threads = int(os.environ.get("THREADS", 1))
    par = int(os.environ.get("PAR", 0))
    ctx.prec = prec; t = time.time(); Q = fr.weil_matrix(K, par, L); tb = time.time() - t
    (l0, e0), (l1, e1) = hp.lowest(Q)
    # lam0_upper: outward-rounded upper end of the Rayleigh-quotient ball (a rigorous upper bound for the sector's lambda0 when
    # the matrix comes from fast_recurrence with rigorous=True, i.e. quadrature error included); lam0_mid is for display only
    r = dict(K=K, par=par, L=L, prec=prec, rigorous_quadrature=True, lam0_upper=(l0.upper() * (1 + arb(10) ** -13)).str(15, radius=False), lam0_mid=l0.mid().str(15, radius=False),
             lam0_rad=l0.rad().str(3, radius=False), res0=e0.upper().str(3, radius=False), lam1=l1.mid().str(8, radius=False), t_build=round(tb), t=round(time.time() - t))
    print(json.dumps(r), flush=True); return r
if __name__ == "__main__":
    out, Ks, Ls, prec = sys.argv[1], [int(k) for k in sys.argv[2].split(",")], sys.argv[3].split(","), int(sys.argv[4])
    jobs = sorted([(K, L, prec) for K in Ks for L in Ls], key=lambda a: -a[0])
    with Pool(min(len(jobs), 16)) as P, open(out, "a") as f:
        for r in P.imap_unordered(job, jobs): f.write(json.dumps(r) + "\n"); f.flush()
