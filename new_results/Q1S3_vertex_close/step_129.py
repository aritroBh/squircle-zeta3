import pickle, mpmath as mp, itertools
from fractions import Fraction as Fr
from q1s_exact import femb, BASISEL
from q1s_symbol import tadd, tscale, t_conj, t_shuffle, t_from_vec, t_tensor, letter_vec
import sympy as sp
mp.mp.dps=80
sols=pickle.load(open('sunit_li3.pkl','rb'))
X=pickle.load(open('q3_Xli2.pkl','rb'))['X']
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def Li2(i): return mp.polylog(2, femb(sols[i]['z'],1))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Xval(V,W):
    v = logu(V)*logu(W)/2
    for i,c in X[(V,W)]: v += mp.mpf(float(c))*Li2(i)
    return v
# artificial V*: 3*(2ReLi3(z145)) + 2*log(a)X(a,p2) - 1*log(q1)X(a,q1)
Vstar = 3*(2*ReLi3(145)) + 2*logu('a')*Xval('a','p2') - logu('q1')*Xval('a','q1')
print("V* =", mp.nstr(Vstar,20))
# its symbol: 3*g145 + 2*(a shuffle (a x p2)) - (q1 shuffle (a x q1))
def dec_sym(u,v,w):
    t={}
    for perm in [(u,v,w),(v,u,w),(v,w,u)]: t[perm]=t.get(perm,Fr(0))+1
    return t
g = tadd(tscale(tadd(sols[145]['S'], t_conj(sols[145]['S'])),Fr(3)),
         tadd(tscale(dec_sym('a','a','p2'),Fr(2)), tscale(dec_sym('q1','a','q1'),Fr(-1))))
# now solve for this symbol in the 55-column basis
B4=('a','p2','q1','q2')
BIDX=list(itertools.product(B4,repeat=3))
sel=pickle.load(open('q3_sel.pkl','rb'))['sel']
from q1s_exact import fconj, feq
pairs=[]; used=set()
for i,s in enumerate(sols):
    if i in used: continue
    j=next((jj for jj,s2 in enumerate(sols) if jj!=i and feq(s2['z'],fconj(s['z']))), None)
    used.add(i)
    if j is not None: used.add(j)
    pairs.append((i,j))
def coord(T):
    v=sp.zeros(64,1)
    for k,c in T.items(): v[BIDX.index(k)]=sp.Rational(c)
    return v
def deco_g(u,v,w):
    t={}
    for perm in [(u,v,w),(v,u,w),(v,w,u)]: t[perm]=t.get(perm,Fr(0))+1
    return coord(t)
Gp={}
for i,j in pairs:
    Gp[(i,j)]=coord(tadd(sols[i]['S'], t_conj(sols[i]['S'])))
cols=[Gp[val] if k=='li3' else deco_g(*val) for k,val in sel]
M=sp.Matrix.hstack(*cols)
t=coord(g)
xsol=M.gauss_jordan_solve(t)[0]
# realize
rec = mp.mpc(0)
for i2,(k,val) in enumerate(sel):
    if xsol[i2]==0: continue
    c=mp.mpf(float(sp.Rational(xsol[i2])))
    if k=='li3':
        rec += c*2*ReLi3(val[0])
    else:
        rec += c*logu(val[0])*Xval(val[1],val[2])
print("recovered:", mp.nstr(rec,20))
print("R* = V* - rec =", mp.nstr(Vstar-rec,10))
