"""Derive log A coefficients from the rational logarithmic derivative, then Bell."""
from pathlib import Path
from math import factorial
import sys,pickle
import sympy as sp
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'TOW3_tower_depth6'))
sys.path.insert(0,str(ROOT/'TOW2_tower_all_depths'))
import tow3lib as T
from genintegrand import S_poly
# S_poly is the finite Stirling-sum definition, independent of the fitted R's.
d,m=T.d,T.m
v,x=sp.symbols('v x'); i=sp.I
G=v**4-6*v**2+1
t=v/4+sp.Rational(1,2)-1/(4*v)
assert t.subs(v,i)==(1+i)/2 and sp.diff(t,v).subs(v,i)==0
z=sp.series(t.subs(v,i+x)-(1+i)/2,x,0,7).removeO()
R1,R2,R3=sp.symbols('R1 R2 R3')
logs=sp.Poly(sp.expand(R1*z+(R2-R1**2)*z**2/2+(R3-3*R1*R2+2*R1**3)*z**3/6),x)
alg=(2*d-1)/v+(2-2*d)/(v+i)-sp.diff(G,v)/G
Ls={}
for r in range(1,7):
    algebra=sp.simplify(sp.diff(alg,v,r-1).subs(v,i))/factorial(r)
    raw=sp.expand(algebra+logs.nth(r))
    print(f'ell_{r} =',sp.expand(raw*factorial(r)),flush=True)
    Ls[r]=sp.cancel(raw.subs({R1:T.T.R[1],R2:T.T.R[2],R3:T.T.R[3]}))
    assert sp.cancel(Ls[r]-T._load_L()[r])==0,('L cache',r)
print('PASS independently derived L1..L6 equal L_parsed.pkl symbolically in d,m',flush=True)
a={0:sp.Integer(1)}
for n in range(1,7):
    a[n]=sp.cancel(sum(r*Ls[r]*a[n-r] for r in range(1,n+1))/n)
    assert sp.cancel(a[n]-T.arho(n))==0,n
print('PASS Bell recurrence a0=1, n*a_n=sum(r*Lr*a_(n-r)): a1..a6 equal submitted arho',flush=True)
# Exact finite-polynomial S derivatives: include q=0, balanced and adjacent keys.
count=0; ts=(1+i)/2
for dd in range(1,10):
    for mm in range((dd+1)//2,dd+1):
        coeffs=S_poly(dd,mm)
        vals=[sp.expand(sum(sp.Rational(c)*sp.ff(n,l)*ts**(n-l) for n,c in enumerate(coeffs) if n>=l)) for l in range(5)]
        for l in range(1,5):
            assert sp.cancel(vals[l]-vals[0]*T.T.R[l].subs({d:dd,m:mm}))==0,(dd,mm,l)
            count+=1
print(f'PASS {count} exact R1..R4 checks from finite S polynomial, d=1..9, every p>=q',flush=True)
# Compare to the exact saved towers: independent initial data ratios.
H=pickle.load(open(ROOT/'TOWN1_tower_large_n'/'data'/'towers_D8_K52.pkl','rb'))
def expr(q): return q[0]+sp.sqrt(2)*q[1]+i*q[2]+i*sp.sqrt(2)*q[3]
counts={5:0,6:0}
for dd in (5,6,7,8,12,20,40,52):
    for mm in sorted({dd,dd-1,(dd+1)//2}):
        key=(mm,dd-mm); base=expr(H[key][0][dd])
        for delta in (5,6):
            expected=a[delta].subs({d:dd,m:mm})
            assert sp.expand(expr(H[key][delta][dd])-base*expected)==0,(key,delta)
            counts[delta]+=1
print('PASS exact initial-tower spot checks',counts,'including d=5, q=0, adjacent and diagonal keys',flush=True)
