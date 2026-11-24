"""Symbolic M3, M5 (2x2) from the BSST associator; compare M3 with SS (3.19)."""
import sys, pickle; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import sympy as sp
from fractions import Fraction as F
import phi
s1,s2,s3,s4,s5=sp.symbols('s1 s2 s3 s4 s5')   # s_i = s_{i,i+1}; s1=s12,s2=s23,s3=s34,s4=s45,s5=s51
def kin_sym(S):
    S1,S2,S3,S4,S5=S
    s13=S4-S1-S2; s24=S5-S2-S3
    return dict(s12=S1,s13=s13,s23=S2,s24=s24,s34=S3,s123=S4,s234=S5)
def sym_cols(S, Wmax):
    k=kin_sym(S); e0,e1=phi.e_mats(**k)
    e0=[[sp.expand(x) for x in r] for r in e0]; e1=[[sp.expand(x) for x in r] for r in e1]
    e={'0':e0,'1':e1}
    cur={'':[sp.Integer(1)]+[sp.Integer(0)]*5}
    res={}
    for W in range(1,Wmax+1):
        nxt={}
        for w,v in cur.items():
            for a in '01':
                M=e[a]; nxt[w+a]=[sp.expand(sum(M[i][j]*v[j] for j in range(6))) for i in range(6)]
        cur=nxt
        if W>=2:
            nb=phi.E.table(W).nullity
            tot=[[0]*nb for _ in range(6)]
            for w,v in cur.items():
                if all(x==0 for x in v): continue
                zv=phi.zeta_vec(w)
                for i in range(6):
                    if v[i]!=0:
                        for kk in range(nb):
                            if zv[kk]!=0: tot[i][kk]+= v[i]*sp.Rational(zv[kk].numerator,zv[kk].denominator)
            res[W]=[[sp.expand(x) for x in r] for r in tot]
    return res
S=(s1,s2,s3,s4,s5)
res=sym_cols(S,5)
def F2_odd(S,k):
    S1,S2=S[0],S[1]
    return -sp.Rational(1,k)*(S1**k+S2**k-(S1+S2)**k)
def swap23(expr):
    # 2<->3 : s12<->s13, s24<->s34, s23,s45,s51 fixed
    return expr.subs({s1:s4-s1-s2, s3:s5-s2-s3},simultaneous=True)
Ms={}
for k in (3,5):
    A=sp.expand(res[k][0][0]+F2_odd(S,k))     # F11 coefficient of zeta_k  (basis idx 0 = zeta_k)
    B=sp.expand(res[k][1][0])                 # F^{(32)} coefficient of zeta_k
    Ms[k]=sp.Matrix([[A,B],[swap23(B),swap23(A)]])
# SS (3.19)
s13=s4-s1-s2; s24=s5-s2-s3
m11=s3*(-s1*(s1+2*s2+s3)+s3*s4+s4**2)+s1*s5*(s1+s5)
m12=-s13*s24*(s1+s2+s3+s4+s5)
m21=s1*s3*(s1+s2+s3-2*(s4+s5))
m22=(s2+s3)*((s1+s2)*(s1+s3)-2*s1*s4)-(2*s1*s3-s24**2+2*s2*(s3+s4))*s5+s4*s5**2
SS3=sp.Matrix([[m11,m12],[m21,m22]])
print('M3 (associator-derived) - M3 (SS 3.19) =',sp.simplify(sp.expand(Ms[3]-SS3)))
print('M3[0,0] =',sp.factor(Ms[3][0,0]))
pickle.dump({k:Ms[k] for k in Ms},open('M35_sym.pkl','wb'))
# also check weight-5 Q: F11's zeta2*zeta3 coefficient etc not needed
