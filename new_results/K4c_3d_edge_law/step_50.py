
import json, mpmath as mp
from mpmath import mpf
mp.mp.dps = 70
from pslq_lib import pslq
from fractions import Fraction

def get_tr(fn):
    J = json.load(open(fn)); S = [mpf(s) for s in J["Sdk"]]
    S0,S1,S2,S3,S4 = S
    L = mpf(-6)/5*mp.log(2) + mpf(13)/5*mp.log(3)
    e1=L; e2=L**2/2-L; e3=L**3/6-L**2+L
    t0=S0; t1=S1-e1*t0; t2=S2-e1*t1-e2*t0; t3=S3-e1*t2-e2*t1-e3*t0
    return t0,t1,t2,t3

pi2=mp.pi**2; z3=mp.zeta(3)
sq2,sq3,sq5,sq6,sq10 = mp.sqrt(2),mp.sqrt(3),mp.sqrt(5),mp.sqrt(6),mp.sqrt(10)

def mk_basis(algos, trans):
    b, n = [], []
    for rn, rv in algos:
        for cn, cv in trans:
            b.append(rv*cv); n.append(f"{rn}*{cn}")
    return b, n

t0,t1,t2,t3 = get_tr("modelA_h96.json")
trans = [("1",mpf(1)),("pi2",pi2),("z3",z3)]
for aname, algos in [("Q(sq2,sq5)", [("1",mpf(1)),("sq2",sq2),("sq5",sq5),("sq10",sq10)]),
                     ("Q(sq2,sq3)", [("1",mpf(1)),("sq2",sq2),("sq3",sq3),("sq6",sq6)]),
                     ("Q(sq2)", [("1",mpf(1)),("sq2",sq2)]),
                     ("Q(sq5)", [("1",mpf(1)),("sq5",sq5)])]:
    B, N = mk_basis(algos, trans)
    for tn, tv in [("t2",t2),("t3",t3)]:
        rel, res = pslq([tv]+B, tol=mpf('1e-46'))
        if rel:
            # verify residual directly
            rr = rel[0]*tv + sum(rel[i+1]*B[i] for i in range(len(B)))
            ok = abs(rr) < mpf('1e-44')
            a0=rel[0]
            combos=[(N[i], Fraction(rel[i+1])/Fraction(-a0)) for i in range(len(B)) if rel[i+1]!=0]
            print(f"{tn} on {aname}: FOUND resid={mp.nstr(res,3)} direct={mp.nstr(rr,3)} {'REAL?' if ok else 'weak'}")
            print("    ", tn, "=", " + ".join(f"({c})*{n}" for n,c in combos))
        else:
            print(f"{tn} on {aname}: none")
