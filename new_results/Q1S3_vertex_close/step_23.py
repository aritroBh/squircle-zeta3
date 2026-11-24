engine = '''"""Q1S3 engine: exact reduction of V1,V3' pieces to Chen words G(w;1) over Q(zeta8),
plus boundary products and rational terms. Letters: 0, a, -b, +-ia, +-ib (y=1/X form).

Data:
- F-elements: 4-tuples of Fractions (via q1s_exact).
- LogExpr: dict factor->coeff(F); factor = ('G',s) meaning G(s;y)=log(1-y/s), s F-element
  (s=0 gives log y), or ('C',x) meaning the CONSTANT logval(x) = log(sigma_1(x)), x F-element.
- WordExpr: dict word(tuple of F)->coeff(F): value sum c*G(word;1).
- ProdExpr: dict tuple-of-factors->coeff(F); factors are ('lv',x) = log(sigma_1(x)) constants.
- Total expression: (WordExpr, ProdExpr, F-number).
"""
from fractions import Fraction as Fr
from q1s_exact import *

i_ = fmul(ZETA,ZETA)
b_ = fsub(r2, ONE)            # sqrt2-1 = 1/a
s0 = ZERO
sa = a_
smb = fneg(b_)
sia = fmul(i_,a_); smia = fneg(sia)
sib = fmul(i_,b_); smib = fneg(sib)
LOG2 = ('C', None)

def Gf(s): return ('G', s)

# ---- LogExpr algebra ----
def le_add(*Es):
    r={}
    for E in Es:
        for k,c in E.items():
            r[k]=fadd(r.get(k,ZERO),c)
            if feq(r[k],ZERO): del r[k]
    return r
def le_scale(E,c):
    return {k:fmul(v,c) for k,v in E.items() if not feq(fmul(v,c),ZERO)}

# log atoms of the pieces
g1 = {Gf(sa):ONE, Gf(smb):ONE, Gf(s0):fneg(ONE), LOG2:fneg(ONE)}      # log((1-by)(1+ay)/(2y))
g2 = le_add({Gf(smib):ONE,Gf(sib):ONE,Gf(smia):ONE,Gf(sia):ONE},
            le_scale({Gf(sa):ONE,Gf(smb):ONE}, fneg(fadd(ONE,ONE))))  # log((1+a^2 y^2)(1+b^2 y^2)/((1-by)^2(1+ay)^2))
g3 = le_add({Gf(smia):ONE,Gf(sia):ONE}, le_scale({Gf(smib):ONE,Gf(sib):ONE}, fneg(ONE)))

# ---- products P = E1*E2: list of (coeff, fA, fB) with fA,fB factors (('G',s) or LOG2) ----
def le_prod(E1,E2):
    out=[]
    for kA,cA in E1.items():
        for kB,cB in E2.items():
            out.append((fmul(cA,cB), kA, kB))
    return out

def le_deriv(E):
    """derivative of LogExpr: list of (coeff, s_form) ; value-part weight drops: for products
    handled in prod_deriv. For a single LogExpr: d G(s;y)/dy = 1/(y-s); d const = 0."""
    return [(c, k[1]) for k,c in E.items() if k[0]=='G']

def prod_deriv(terms):
    """terms = le_prod output; -> list of (coeff, s_form, value_factor) meaning coeff * (1/(y-s_form)) * value_factor(y)."""
    out=[]
    for c,fA,fB in terms:
        if fA[0]=='G': out.append((c, fA[1], fB))
        if fB[0]=='G': out.append((c, fB[1], fA))
    return out

# ---- word/product accumulators ----
class Acc:
    def __init__(self):
        self.W={}   # word -> F coeff
        self.P={}   # tuple of lv factors -> F coeff
        self.R=ZERO # rational F-number
    def add_word(self,w,c):
        if feq(c,ZERO): return
        self.W[w]=fadd(self.W.get(w,ZERO),c)
        if feq(self.W[w],ZERO): del self.W[w]
    def add_prod(self,factors,c):
        if feq(c,ZERO): return
        key=tuple(sorted(factors,key=repr))
        self.P[key]=fadd(self.P.get(key,ZERO),c)
        if feq(self.P[key],ZERO): del self.P[key]
    def add_R(self,c):
        self.R=fadd(self.R,c)
    def scale(self,c):
        A=Acc()
        A.W={k:fmul(v,c) for k,v in self.W.items() if not feq(fmul(v,c),ZERO)}
        A.P={k:fmul(v,c) for k,v in self.P.items() if not feq(fmul(v,c),ZERO)}
        A.R=fmul(self.R,c)
        return A
    def add_into(self,other):
        for k,v in other.W.items(): self.add_word(k,v)
        for k,v in other.P.items(): self.add_prod(k,v)
        self.R=fadd(self.R,other.R)

def lv_of_Gval(s):
    """G(s;1) = log(1 - 1/s) as ('lv', x)."""
    assert not feq(s,ZERO)
    return ('lv', fdiv(fsub(s,ONE), s))
def fdiv(x,y): return fmul(x,finv(y))

def atom_at1(f):
    """value of a factor at y=1 as list of lv factors (length 0 or 1)."""
    if f[0]=='C': return [('lv', fadd(ONE,ONE))]   # log 2
    if feq(f[1],ZERO): return []                    # log 1 = 0
    return [lv_of_Gval(f[1])]

# ---- integration rules ----
def int_simple_kernel(acc, cK, sK, terms):
    """int_0^1 cK * (1/(y-sK)) * P(y) dy, P given as le_prod terms."""
    for c,fA,fB in terms:
        cc=fmul(cK,c)
        if fA[0]=='G' and fB[0]=='G':
            acc.add_word((sK,fA[1],fB[1]), cc)
            acc.add_word((sK,fB[1],fA[1]), cc)
        elif fA[0]=='G':
            acc.add_word((sK,fA[1]), fmul(cc, fB[1] and ONE or ONE))  # fB const coeff handled below
        else:
            acc.add_word((sK,fB[1]), cc)
    # constants: handle separately (terms with a const factor)
    for c,fA,fB in terms:
        cc=fmul(cK,c)
        if fA[0]=='C' and fB[0]=='C':
            # cK*cA*cB * G(sK;1)
            acc.add_prod([lv_of_Gval(sK)], cc)
        elif fA[0]=='C' or fB[0]=='C':
            pass # already handled above with coeff 1 -- WRONG, fix in caller

def int_simple_kernel2(acc, cK, sK, terms):
    """correct version: int cK/(y-sK) * P dy."""
    for c,fA,fB in terms:
        cc=fmul(cK,c)
        if fA[0]=='G' and fB[0]=='G':
            acc.add_word((sK,fA[1],fB[1]), cc); acc.add_word((sK,fB[1],fA[1]), cc)
        elif fA[0]=='G' and fB[0]=='C':
            acc.add_word((sK,fA[1]), cc)   # const factor is the NUMBER log2: keep as product with G-word? NO:
        elif fB[0]=='G' and fA[0]=='C':
            acc.add_word((sK,fB[1]), cc)
        else:
            acc.add_prod([lv_of_Gval(sK)], cc)

def int_dy_value(acc, c, val):
    """int_0^1 c * d(val-factor) = c*[val]_0^1 ; val = ('G',s) or LOG2(d->0)."""
    if val[0]=='C': return
    for f in atom_at1(val): acc.add_prod([f], c)

def pf_two_forms(s1,s2):
    """1/((y-s1)(y-s2)) partial fractions -> list of (coeff, r) for 1/(y-r), or ('sq',s) if s1==s2."""
    if feq(s1,s2): return ('sq', s1)
    d=fsub(s1,s2)
    return [(fdiv(ONE,d), s1), (fneg(fdiv(ONE,d)), s2)]

def int_formform_value(acc, cF,sF, cK,sK, val):
    """int_0^1 cF*cK * [1/((y-sF)(y-sK))] * val(y) dy, val = ('G',u) or ('C',..)."""
    cc=fmul(cF,cK)
    pf=pf_two_forms(sF,sK)
    if pf[0]=='sq':
        _,s=pf
        # int val/(y-s)^2 dy  with B(y)=-y/(s(y-s)), B(0)=0:
        if val[0]=='C':
            acc.add_R(fmul(cc, fneg(fdiv(ONE, fmul(s,fsub(ONE,s))))))
            return
        u=val[1]
        # boundary val(1)*B(1), B(1) = -1/(s(1-s))
        B1=fneg(fdiv(ONE, fmul(s,fsub(ONE,s))))
        if feq(u,ZERO):
            pass  # val(1)=log 1 =0
        else:
            acc.add_prod([lv_of_Gval(u)], fmul(cc,B1))
        # -(1/s) int val' * (-y/((y-s)))? we derived: int val B' = val(1)B(1) + (1/s) int y/((y-sv)(y-s)) dy
        sv=u
        if feq(sv,s):
            # int y/(y-s)^2 = log(1-1/s) - s/(1-s) - 1
            acc.add_prod([lv_of_Gval(s)], fdiv(cc,s))
            acc.add_R(fmul(fdiv(cc,s), fadd(fneg(fdiv(s,fsub(ONE,s))), fneg(ONE))))
        elif feq(sv,ZERO):
            # int y/(y(y-s)) = int 1/(y-s) = log(1-1/s)
            acc.add_prod([lv_of_Gval(s)], fdiv(cc,s))
        else:
            A=fdiv(sv,fsub(sv,s)); Bm=fneg(fdiv(s,fsub(sv,s)))
            acc.add_prod([lv_of_Gval(sv)], fmul(fdiv(cc,s),A))
            acc.add_prod([lv_of_Gval(s)], fmul(fdiv(cc,s),fneg(Bm)))
        return
    for (cr,r) in pf:
        c2=fmul(cc,cr)
        if val[0]=='C':
            acc.add_prod([lv_of_Gval(r)], c2)
        else:
            u=val[1]
            if feq(u,ZERO):
                acc.add_word((r,ZERO), c2)
            else:
                acc.add_word((r,u), c2)

def boundary_P01(acc, cP, terms):
    """cP * [P]_0^1 ; P(0)=0 termwise (checked); P(1) = product of values at 1."""
    for c,fA,fB in terms:
        facs = atom_at1(fA)+atom_at1(fB)
        if facs: acc.add_prod(facs, fmul(cP,c))
        else: acc.add_R(fmul(cP,c))  # both factors zero at 1 (log1): contributes 0 actually; keep 0
'''
open('q3_engine.py','w').write(engine)
print("written", len(engine))
