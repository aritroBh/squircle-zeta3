import pickle, mpmath as mp
from q1s_exact import femb
from q1s_symbol import tadd, tscale, t_shuffle, t_tensor, t_from_vec, letter_vec
from fractions import Fraction as Fr
import sympy as sp, itertools
mp.mp.dps=60
BAS={'a':None,'p2':None,'q1':None,'q2':None}
from q1s_exact import BASISEL
def logu(u): return mp.log(femb(BASISEL[u],1))
def Gvw(v,w):
    vn=femb(BASISEL[v],1); wn=femb(BASISEL[w],1)
    return mp.quad(lambda t: mp.log(1-t/wn)/(t-vn), [0,1])
# convention check: S(log u * G(v,w;1)) =? u shuffle (v x w)
u,v,w = 'a','p2','q1'
LHS_sym = t_shuffle(t_from_vec(letter_vec(BASISEL[u])), t_tensor(t_from_vec(letter_vec(BASISEL[v])), t_from_vec(letter_vec(BASISEL[w]))))
# generator built in deco_gens:
t={}
for perm in [(u,v,w),(v,u,w),(v,w,u)]: t[perm]=t.get(perm,Fr(0))+1
print("symbol match:", LHS_sym==t)
sol = pickle.load(open('q3_solve.pkl','rb'))
# decomp value of V3p
tot=mp.mpc(0)
for (u,v,w),c in sol['V3p']['dec']:
    val = logu(u)*Gvw(v,w)
    tot += mp.mpf(float(c))*val
    print("  (%s,%s,%s) c=%-4s log(u)G=%s"%(u,v,w,c,mp.nstr(val,12)))
print("decomp total:", mp.nstr(tot,16))
