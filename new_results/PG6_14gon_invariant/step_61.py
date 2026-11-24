import pickle, sys
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import engine14 as E
# FT element = (A,B,C,d): (A+B a+C a^2)/d. Letters: 2, alpha, u=alpha^2+alpha-1.
# factor: z = sign * 2^e0 * alpha^e1 * u^e2. Everything has norm +/-2^k.
# Norm of alpha: minpoly x^3+x^2-2x-1, product of roots = -(-1)=1 -> N(alpha)= -1? for monic cubic x^3+px^2+qx+r, N = -r = 1. Actually N = (-1)^3 * r = -r. r=-1 -> N=1? 
# Let me just compute norms via resultant.
import sympy as sp
def ft_norm(x):
    A,B,C,d = x
    t=sp.symbols('t')
    f=sp.Rational(A)+sp.Rational(B)*t+sp.Rational(C)*t**2
    g=t**3+t**2-2*t-1
    return Fr(sp.resultant(f,g,t))/(d**3)
AL=(Fr(0),Fr(1),Fr(0),1)
UU=(Fr(-1),Fr(1),Fr(1),1)  # u = alpha^2+alpha-1
print("N(2)=8, N(alpha)=", ft_norm(AL), " N(u)=", ft_norm(UU))
# So N(alpha), N(u) are units (+/-1). N(2)=8. exponent of 2 in z from norm.
import math
LOG2=math.log(2); EMB=[2*math.cos(2*math.pi/7),2*math.cos(4*math.pi/7),2*math.cos(6*math.pi/7)]
UEMB=[e*e+e-1 for e in EMB]
def ft_embeds(x):
    A,B,C,d=x
    A,B,C=float(A),float(B),float(C)
    return [(A+B*e+C*e*e)/d for e in EMB]
def factor_ft(x, bound=20):
    """find (sign,e0,e1,e2) with z=sign*2^e0*alpha^e1*u^e2, or None."""
    if x[0]==0 and x[1]==0 and x[2]==0: return None
    N=ft_norm(x)
    # e0 from norm: N = +/- 2^(3 e0) (since N(alpha)=N(u)=+/-1, N(2)=8)
    an=abs(Fr(N))
    e=0; v=Fr(1)
    # count power of 2
    num,den=an.numerator, an.denominator
    e0n=0
    while num%2==0 and num!=0: num//=2; e0n+=1
    e0d=0
    while den%2==0 and den!=0: den//=2; e0d+=1
    if num!=1 or den!=1: return ('NOT2', num, den)
    e3 = e0n-e0d  # = 3*e0
    if e3%3!=0: return ('BAD3', e3)
    e0 = e3//3
    # now z/2^e0 = +/- alpha^e1 u^e2 (a unit). Determine e1,e2 via embeddings (logs of abs).
    zv = ft_embeds(x)
    zv = [zz/(2.0**e0) for zz in zv]
    la=[math.log(abs(EMB[i])) for i in range(3)]
    lu=[math.log(abs(UEMB[i])) for i in range(3)]
    lz=[math.log(abs(zz)) for zz in zv]
    # solve lz = e1*la + e2*lu (3 eqns, 2 unknowns) least squares
    import numpy as np
    M=np.array([[la[i],lu[i]] for i in range(3)]); b=np.array(lz)
    sol,res,rank,sv=np.linalg.lstsq(M,b,rcond=None)
    e1,e2=int(round(sol[0])),int(round(sol[1]))
    # verify exactly
    cand = E.fpow(AL,e1) if e1>=0 else E.fpow(E.finv(AL),-e1)
    cand = E.fmul(cand, E.fpow(UU,e2) if e2>=0 else E.fpow(E.finv(UU),-e2))
    cand = E.fs(cand, Fr(2)**e0)
    ratio = E.fdiv(x, cand)
    if E.feq0(E.fsub(ratio,E.ONE)): return (1,e0,e1,e2)
    if E.feq0(E.fadd(ratio,E.ONE)): return (-1,e0,e1,e2)
    return ('UNITFAIL',(e0,e1,e2),ratio)
# test against letters14 facs
L = pickle.load(open('letters14.pkl','rb'))
g = pickle.load(open('geom14.pkl','rb'))
# build some quantities: x_i, 1-x_i
tests={}
for i in [1,2,5,6]:
    tests['x%d'%i]=g['X'][i]
    tests['1-x%d'%i]=E.fsub(E.ONE,g['X'][i])
    tests['1/x%d'%i]=E.finv(g['X'][i])
    tests['(x%d-1)/x%d'%(i,i)]=E.fdiv(E.fsub(g['X'][i],E.ONE),g['X'][i])
for nm,z in tests.items():
    print(nm, factor_ft(z))
