import sys
sys.path.insert(0,'.')
from oct9c_common import *
D = load_all()
mp.mp.dps = 100
e = 1
eng = get_engine(e)
W3 = D['s123'][e]['W3']
NUM = atom_values(eng, e, D['lam'], 100)
W3val = mp.mpf(sp.N(W3.subs(list(NUM.items())), 100))
Gargs, Dargs = rhs_terms_235()
sol = D['s123'][e]['sol']
RHS = eval_rhs_235(sol, Gargs, Dargs, 100)
R = W3val - RHS
print("W3val =", mp.nstr(W3val, 30))
print("RHS   =", mp.nstr(RHS, 30))
print("R     =", mp.nstr(R, 30))
z2 = mp.pi**2/6
inv = [mp.zeta(3), z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)]
rels = mp.pslq([R]+inv, maxcoeff=10**12, maxsteps=10000)
print("pslq:", rels)
if rels:
    err = abs(rels[0]*R + sum(r_*v for r_,v in zip(rels[1:], inv)))
    print("err:", mp.nstr(err,3))
