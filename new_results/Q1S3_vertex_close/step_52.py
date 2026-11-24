import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp
mp.mp.dps=80
r2=mp.sqrt(2); a=1+r2; b=r2-1
ge=GEval(80, M=48, gln=5)
def acc_value(acc):
    tot=mp.mpc(0)
    for (w,lvm),c in acc.W.items():
        wn=tuple(femb(s,1) for s in w)
        v=ge.word(wn)
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    for lvm,c in acc.P.items():
        v=mp.mpc(1)
        for f in lvm: v*=mp.log(femb(f[1],1))
        tot+=femb(c,1)*v
    tot+=femb(acc.R,1)
    return tot
# single term: int (1/y) log(1-by) log(1+iby) dy ; log(1+iby) has letter -ia
P = le_prod({Gf(sa):ONE},{Gf(smia):ONE})
print("P terms:", P)
acc=Acc(); int_simple_kernel(acc, ONE, s0, P)
vd = mp.quad(lambda y: mp.log(1-b*y)*mp.log(1+1j*b*y)/y, [0,1])
print("single term: engine", mp.nstr(acc_value(acc),18), " direct", mp.nstr(vd,18))
# with log y factor: int (1/y) log y log(1-by) dy
P2 = le_prod({Gf(s0):ONE},{Gf(sa):ONE})
acc2=Acc(); int_simple_kernel(acc2, ONE, s0, P2)
vd2 = mp.quad(lambda y: mp.log(y)*mp.log(1-b*y)/y, [0,1])
print("logy term: engine", mp.nstr(acc_value(acc2),18), " direct", mp.nstr(vd2,18))
# with const: int (1/y) log2 * log(1-by) dy
P3 = le_prod({LOG2:ONE},{Gf(sa):ONE})
acc3=Acc(); int_simple_kernel(acc3, ONE, s0, P3)
vd3 = mp.quad(lambda y: mp.log(2)*mp.log(1-b*y)/y, [0,1])
print("const term: engine", mp.nstr(acc_value(acc3),18), " direct", mp.nstr(vd3,18))
