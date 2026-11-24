"""Q1S symbol engine: symbols of iterated integrals in V^{otimes n}, V = (O_{F,S}^x/tors)(x)Q, F=Q(zeta8).
Tensor: dict n-tuple of basis letters ('a','p2','q1','q2') -> Fraction.
Vector: dict basis letter -> Fraction.
Symbol via Goncharov coproduct recursion with shuffle-regularized weight-1 letters.
"""
from fractions import Fraction as Fr
from itertools import combinations
from q1s_exact import *

B=('a','p2','q1','q2')
def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def vsub(a,b): return vadd(a,{k:-c for k,c in b.items()})
def vmul(v,a):
    a=Fr(a); return {k:c*a for k,c in v.items() if c*a!=0}

def tadd(*Ts):
    r={}
    for T in Ts:
        for k,c in T.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
def tscale(T,a):
    a=Fr(a); return {k:c*a for k,c in T.items() if c*a!=0}
def t_from_vec(v):
    return {(k,):c for k,c in v.items()}
def t_tensor(T1,T2):
    r={}
    for k1,c1 in T1.items():
        for k2,c2 in T2.items():
            r[k1+k2]=r.get(k1+k2,Fr(0))+c1*c2
    return {k:c for k,c in r.items() if c}
def t_shuffle(T1,T2):
    # shuffle product of tensors (multilinear): sum over interleavings
    r={}
    for k1,c1 in T1.items():
        for k2,c2 in T2.items():
            n,m=len(k1),len(k2)
            for pos in combinations(range(n+m),n):
                p1=set(pos); w=[]; i=j=0
                for t in range(n+m):
                    if t in p1: w.append(k1[i]); i+=1
                    else: w.append(k2[j]); j+=1
                w=tuple(w)
                r[w]=r.get(w,Fr(0))+c1*c2
    return {k:c for k,c in r.items() if c}

def letter_vec(elem):
    f=factor_sunit(elem)
    assert f is not None, ("not an S-unit", fshow(elem))
    return f[1]

def sym_w1(a0,a1,a2):
    """regularized symbol of I(a0;a1;a2)=log((a2-a1)/(a0-a1)): drop zero factors."""
    v={}
    num=fsub(a2,a1); den=fsub(a0,a1)
    if not feq(num,ZERO): v=vadd(v,letter_vec(num))
    if not feq(den,ZERO): v=vsub(v,letter_vec(den))
    return v

_MEMO={}
def sym_I(word, a0, an1):
    """symbol of I(a0; word; an1), word = tuple of F-elements (letters).
    Rightmost-iteration recursion: S(I_n) = sum_j [I(a_{j-1};a_j;a_{j+1})] x S(I minus letter j).
    (Validated on controls: Li3, Li2, Landen, Li21, shuffles, integrability.)"""
    n=len(word)
    key=(word,a0,an1)
    if key in _MEMO: return _MEMO[key]
    if n==1:
        r=t_from_vec(sym_w1(a0,word[0],an1)); _MEMO[key]=r; return r
    pts=(a0,)+tuple(word)+(an1,)
    res={}
    for j in range(1,n+1):
        lv=sym_w1(pts[j-1], pts[j], pts[j+1])
        sub=word[:j-1]+word[j:]
        res=tadd(res, t_tensor(t_from_vec(lv), sym_I(sub, a0, an1)))
    _MEMO[key]=res
    return res

def t_conj(T):
    """complex conjugation on tensor factors: q1<->q2."""
    r={}
    sw={'a':'a','p2':'p2','q1':'q2','q2':'q1'}
    for k,c in T.items():
        r[tuple(sw[x] for x in k)]=r.get(tuple(sw[x] for x in k),Fr(0))+c
    return {k:c for k,c in r.items() if c}
def teq(A,B):
    keys=set(A)|set(B)
    return all(A.get(k,Fr(0))==B.get(k,Fr(0)) for k in keys)
def tsub(A,B): return tadd(A,tscale(B,-1))
