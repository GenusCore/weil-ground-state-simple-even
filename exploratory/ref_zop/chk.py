import os, sys; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "code"))  # for convergence.py
import mpmath as mp, numpy as np
from convergence import galerkin
mp.mp.dps = 120
L, D = float(sys.argv[1]), int(sys.argv[2])
Q = galerkin(L, D)
A = mp.matrix(D, D)
for i in range(D):
    for j in range(D): A[i, j] = mp.mpf(Q[i, j].mid().str(110, radius=False))
def e(m, x): return mp.sqrt((2*m+1)/(2*mp.mpf(L))) * mp.legendre(m, x/L)
cc = [mp.quad(lambda x: e(m, x)*mp.cosh(x/2), [-L, 0, L]) for m in range(D)]
ss = [mp.quad(lambda x: e(m, x)*mp.sinh(x/2), [-L, 0, L]) for m in range(D)]
def sub(Mat, par):
    idx = list(range(par, D, 2)); n=len(idx); S=mp.matrix(n,n)
    for a in range(n):
        for b in range(n): S[a,b]=Mat[idx[a],idx[b]]
    return idx, S
ie, Qe = sub(A,0); io, Qo = sub(A,1)
ce = mp.matrix([cc[i] for i in ie]); so = mp.matrix([ss[i] for i in io])
Me = Qe - 2*ce*ce.T; Mo = Qo + 2*so*so.T
mu_e, Ve = mp.eigsy(Me); mu_o, Vo = mp.eigsy(Mo)
lq_e,_ = mp.eigsy(Qe); lq_o,_ = mp.eigsy(Qo)
a = [ (Ve[:,k].T*ce)[0] for k in range(4)]
b = [ (Vo[:,k].T*so)[0] for k in range(4)]
print("L",L,"D",D)
print("mu_e", [mp.nstr(mu_e[k],4) for k in range(4)], "a_k", [mp.nstr(x,4) for x in a])
print("mu_o", [mp.nstr(mu_o[k],4) for k in range(4)], "b_k", [mp.nstr(x,4) for x in b])
print("lam_e", [mp.nstr(lq_e[k],4) for k in range(3)], "lam_o", [mp.nstr(lq_o[k],4) for k in range(3)])
l0 = lq_e[0]
sec = 2*sum(b[k]**2/(mu_o[k]-l0) for k in range(len(io)) for b_ in [None]) if False else None
bb = [ (Vo[:,k].T*so)[0] for k in range(len(io))]
sec = 2*sum(bb[k]**2/(mu_o[k]-l0) for k in range(len(io)))
print("odd secular 2<s,(Mo-l0e)^-1 s> =", mp.nstr(sec,8), " (need <1)")
# relative weight of the mu_o[0] term
print("term0 share", mp.nstr(2*bb[0]**2/(mu_o[0]-l0)/sec,6))
print("|s|^2*2=", mp.nstr(2*(so.T*so)[0],6))
