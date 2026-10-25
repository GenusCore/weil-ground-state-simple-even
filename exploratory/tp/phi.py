# Rank-one intertwining test.  Q_o d - d Q_e = f(L) w on even f  ==>  nonzero odd eigenvalues = zeros of
#   Phi(z) = sum_k c_k/(a_k - z),  c_k = e_k(L) * int e_k  (even eigenpairs (a_k,e_k)),  sum c_k = 1.
# If all c_k > 0: a_0 < o_0 < a_1 < o_1 < ...  (simple, even ground state + parity alternation).
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))  # for convergence.py
import mpmath as mp
from convergence import galerkin
mp.mp.dps = 140
L, D, nshow = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
Q = galerkin(L, D)
def sector(idx):
    A = mp.matrix([[mp.mpf(Q[i, j].mid().str(130, radius=False)) for j in idx] for i in idx])
    e, V = mp.eigsy(A); o = sorted(range(len(idx)), key=lambda k: e[k])
    return [e[k] for k in o], [[V[i, k] for i in range(len(idx))] for k in o]
ae, Ve = sector(list(range(0, D, 2))); ao, _ = sector(list(range(1, D, 2)))
# e_m(L) = sqrt((2m+1)/2L) P_m(1) = sqrt((2m+1)/2L);  int e_m = sqrt(2L) delta_m0
c = [sum(v[i] * mp.sqrt((4 * i + 1) / (2 * L)) for i in range(len(v))) * v[0] * mp.sqrt(2 * L) for v in Ve]
print(f"L={L} D={D} sum c_k={mp.nstr(sum(c),8)}")
for k in range(nshow):
    print(f" k={k} a_k={mp.nstr(ae[k],6)} c_k={mp.nstr(c[k],5)}  o_k={mp.nstr(ao[k],6)} Phi(o_k)/max|term|={mp.nstr(sum(ci/(a-ao[k]) for ci,a in zip(c,ae))/max(abs(ci/(a-ao[k])) for ci,a in zip(c,ae)),3)}")
