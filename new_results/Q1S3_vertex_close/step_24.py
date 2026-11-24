engine = '''"""Q1S3 engine: exact reduction of the V1,V3' pieces to Chen words G(w;1) over Q(zeta8)
(y=1/X form; letters 0,a,-b,+-ia,+-ib), plus boundary products of log-constants,
plus rational F-terms.

Conventions:
- F-elements: 4-tuples of Fractions (q1s_exact).
- LogExpr: dict factor->F-coeff; factor = ('G',s): G(s;y)=log(1-y/s) (s=0 -> log y);
  ('C',0): the constant log 2.
- Expression (class Acc): 
    W: dict (word, lvmonomial) -> F-coeff : value = c * G(word;1) * prod lvmonomial
    P: dict lvmonomial -> F-coeff          : products of log-constants
    R: F-element                           : rational term
  lvmonomial = sorted tuple of ('lv',x), meaning log(sigma_1(x)) (principal complex log).
"""
from fractions import Fraction as Fr
from q1s_exact import *

i_ = fmul(ZETA,ZETA)
b_ = fsub(r2, ONE)
s0 = ZERO
sa = a_
smb = fneg(b_)
sia = fmul(i_,a_); smia = fneg(sia)
sib = fmul(i_,b_); smib = fneg(sib)
LOG2 = ('C',0)
TWO = fadd(ONE,ONE)

def Gf(s): return ('G', s)
def fdiv(x,y): return fmul(x,finv(x) if False else fmul(x,finv(y)))
def fdiv(x,y): return fmul(x,finv(y))

def le_add(*Es):
    r={}
    for E in Es:
        for k,c in E.items():
            r[k]=fadd(r.get(k,ZERO),c)
            if feq(r[k],ZERO): del r[k]
    return r
def le_scale(E,c):
    return {k:fmul(v,c) for k,v in E.items() if not feq(fmul(v,c),ZERO)}

# log expressions of the pieces (exact algebra; certified later)
g1 = {Gf(sa):ONE, Gf(smb):ONE, Gf(s0):fneg(ONE), LOG2:fneg(ONE)}
g2 = le_add({Gf(smib):ONE,Gf(sib):ONE,Gf(smia):ONE,Gf(sia):ONE},
            le_scale({Gf(sa):ONE,Gf(smb):ONE}, fneg(TWO)))
g3 = le_add({Gf(smia):ONE,Gf(sia):ONE}, le_scale({Gf(smib):ONE,Gf(sib):ONE}, fneg(ONE)))

def le_prod(E1,E2):
    return [(fmul(cA,cB), kA, kB) for kA,cA in E1.items() for kB,cB in E2.items()]

def le_deriv_terms(E):
    """d/dy of LogExpr as list of (coeff, s_form): only 'G' factors contribute."""
    return [(c, k[1]) for k,c in E.items() if k[0]=='G']

def prod_deriv(terms):
    """derivative of le_prod(E1,E2): list of (coeff, s_form, val_factor)."""
    out=[]
    for c,fA,fB in terms:
        if fA[0]=='G': out.append((c, fA[1], fB))
        if fB[0]=='G': out.append((c, fB[1], fA))
    return out

class Acc:
    def __init__(self):
        self.W={}; self.P={}; self.R=ZERO
    def add_word(self,w,lvm,c):
        if feq(c,ZERO): return
        key=(w,lvm)
        self.W[key]=fadd(self.W.get(key,ZERO),c)
        if feq(self.W[key],ZERO): del self.W[key]
    def add_prod(self,lvm,c):
        if feq(c,ZERO): return
        lvm=tuple(sorted(lvm,key=repr))
        self.P[lvm]=fadd(self.P.get(lvm,ZERO),c)
        if feq(self.P[lvm],ZERO): del self.P[lvm]
    def add_R(self,c):
        self.R=fadd(self.R,c)
    def scale(self,c):
        A=Acc()
        for k,v in self.W.items(): A.W[k]=fmul(v,c)
        for k,v in self.P.items(): A.P[k]=fmul(v,c)
        A.R=fmul(self.R,c)
        return A
    def add_into(self,o):
        for k,v in o.W.items():
            self.W[k]=fadd(self.W.get(k,ZERO),v)
            if feq(self.W[k],ZERO): del self.W[k]
        for k,v in o.P.items(): self.add_prod(k,v)
        self.R=fadd(self.R,o.R)

def lv_G1(s):
    """G(s;1) = log(1-1/s) as lv factor."""
    assert not feq(s,ZERO)
    return ('lv', fdiv(fsub(s,ONE), s))
def atom_at1(f):
    """value of factor f at y=1: list of lv factors ([] if the value is 0)."""
    if f[0]=='C': return [('lv',TWO)]
    if feq(f[1],ZERO): return []
    return [lv_G1(f[1])]

# ---------- rules ----------
def int_simple_kernel(acc, cK, sK, terms):
    """int_0^1 cK/(y-sK) * P(y) dy."""
    for c,fA,fB in terms:
        cc=fmul(cK,c)
        if fA[0]=='G' and fB[0]=='G':
            acc.add_word((sK,fA[1],fB[1]), (), cc)
            acc.add_word((sK,fB[1],fA[1]), (), cc)
        elif fA[0]=='G':
            acc.add_word((sK,fA[1]), tuple(atom_at1(fB)), cc)
        elif fB[0]=='G':
            acc.add_word((sK,fB[1]), tuple(atom_at1(fA)), cc)
        else:
            acc.add_prod(tuple(atom_at1(fA)+[lv_G1(sK)]), cc)

def int_dy(acc, c, val):
    """int_0^1 c * d(val) : val a single factor."""
    if val[0]=='C': return
    if feq(val[1],ZERO): acc.add_R(c)   # int c/y dy -- DIVERGENT; must never happen
    else: acc.add_prod(tuple(atom_at1(val)), c)

def pf2(s1,s2):
    if feq(s1,s2): return ('sq',s1)
    d=fsub(s1,s2)
    return [(fdiv(ONE,d),s1),(fneg(fdiv(ONE,d)),s2)]

def int_formform_val(acc, cF,sF, cK,sK, val):
    """int_0^1 cF*cK/((y-sF)(y-sK)) * val(y) dy, val=('G',u) or ('C',0)."""
    cc=fmul(cF,cK)
    pf=pf2(sF,sK)
    if pf[0]=='sq':
        s=pf[1]
        if val[0]=='C':
            acc.add_R(fmul(cc, fneg(fdiv(ONE, fmul(s,fsub(ONE,s)))))); return
        u=val[1]
        B1=fneg(fdiv(ONE, fmul(s,fsub(ONE,s))))
        for f in atom_at1(val): acc.add_prod((f,), fmul(cc,B1))
        if feq(u,s):
            acc.add_prod((lv_G1(s),), fdiv(cc,s))
            acc.add_R(fmul(fdiv(cc,s), fneg(fadd(fdiv(s,fsub(ONE,s)),ONE))))
        elif feq(u,ZERO):
            acc.add_prod((lv_G1(s),), fdiv(cc,s))
        else:
            A=fdiv(u,fsub(u,s)); Bm=fneg(fdiv(s,fsub(u,s)))
            acc.add_prod((lv_G1(u),), fmul(fdiv(cc,s),A))
            acc.add_prod((lv_G1(s),), fmul(fdiv(cc,s),fneg(Bm)))
        return
    for cr,r in pf:
        c2=fmul(cc,cr)
        if val[0]=='C':
            acc.add_prod((lv_G1(r),), c2)
        else:
            u=val[1]
            acc.add_word((r,u), (), c2)

def int_formP(acc, cK, sK, dterms):
    """int_0^1 cK/(y-sK) * P'(y) dy, dterms from prod_deriv: (c,s_form,val)."""
    for c,sF,val in dterms:
        int_formform_val(acc, c,sF, cK,sK, val)

def int_P_dy(acc, c, dterms):
    """int_0^1 c * P'(y) dy = c*[P]_0^1, computed termwise: each (c',s_form,val): int c c' val dG(s_form)."""
    for c2,sF,val in dterms:
        cc=fmul(c,c2)
        if val[0]=='C':
            int_dy(acc, cc, ('G',sF))
        else:
            # int cc * val * dG(sF) : if val=('G',u): word (sF,u) [outer form sF]
            u=val[1]
            acc.add_word((sF,u), (), cc)

def boundary_P_at1(acc, cP, terms):
    """cP * P(1): P(1) = sum of products of factor-values at 1."""
    for c,fA,fB in terms:
        facs=tuple(atom_at1(fA)+atom_at1(fB))
        if facs: acc.add_prod(facs, fmul(cP,c))

# y/((y-r)(y-s)) integrals (for H4,H5): returns Acc with products/rational
def int_y_rr(acc, cc, r, s):
    """int_0^1 cc * y/((y-r)(y-s)) dy."""
    if feq(r,s):
        # y/(y-r)^2 = 1/(y-r) + r/(y-r)^2 : int = log(1-1/r) + r*(-1/(1-r) - 1/r)
        acc.add_prod((lv_G1(r),), cc)
        acc.add_R(fmul(cc, fsub(fneg(fdiv(r,fsub(ONE,r))), ONE)))
    elif feq(r,ZERO):
        acc.add_prod((lv_G1(s),), cc)
    elif feq(s,ZERO):
        acc.add_prod((lv_G1(r),), cc)
    else:
        A=fdiv(r,fsub(r,s)); Bm=fneg(fdiv(s,fsub(r,s)))
        acc.add_prod((lv_G1(r),), fmul(cc,A))
        acc.add_prod((lv_G1(s),), fmul(cc,Bm))
'''
open('q3_engine.py','w').write(engine)
print("ok", len(engine))
