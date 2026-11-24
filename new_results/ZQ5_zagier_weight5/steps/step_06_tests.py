
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import itertools, math, time
import mpmath as mp
import sympy as sp
import w5lib as W
from w5lib import QD

print("=== QD identities ===")
w2=QD(0,1,2); b=w2-1
print("(1+w2)^2 == 3+2w2:", (1+w2)**2 == QD(3,2,2))
phi=QD(Fraction(1,2),Fraction(1,2),5)
print("phi^2 == phi+1:", phi**2 == phi+1)
r3=QD(0,1,3)
print("(1+r3)^2 == 2*(2+r3):", (1+r3)**2 == 2*(2+r3))
print("b^4 == 17-12w2:", b**4 == QD(17,-12,2))
print("sqrt5-2 == phi^-3:", QD(-2,1,5) == phi**-3)
print("3-sqrt5 == 2*phi^-2:", QD(3,-1,5) == 2*(phi**-2))

print("=== factor_sunit: 20 random S-units per D, reconstruct ===")
import random
random.seed(12345)
for D in (8,5,12,13):
    C=W.letters_info(D); d=len(C['letters']); ok=0
    for trial in range(20):
        exps=[random.randint(-6,6) for _ in range(d)]
        if all(e==0 for e in exps): exps[0]=1
        x=QD(1,0,C['letters'][0].D)
        for g,e in zip(C['letters'],exps): x=x*(g**e)
        if random.choice([0,1]): x=-x
        f=W.factor_sunit(x,C)
        assert f==exps, (D,exps,f)
        # reconstruct |x| from coords
        y=QD(1,0,C['letters'][0].D)
        for g,e in zip(C['letters'],f): y=y*(g**e)
        assert y==x or y==-x
        ok+=1
    # rejection tests
    if D==8: bad=[QD(3,1,2)]            # N=7
    if D==5: bad=[QD(3,2,5)]            # N=-11
    if D==12: bad=[QD(4,1,3)]          # N=13
    if D==13: bad=[QD(1,1,13)]          # N=-12= -2^2*3 -> 3 not allowed
    for x in bad:
        assert W.factor_sunit(x,C) is None, (D,x)
    print(f"D={D}: 20/20 random S-units factored+reconstructed, rejections ok; {C['Sdesc']}")

print("=== saturation & kernel on planted matrices ===")
# planted: A (4x2) rank2, left kernel spanned by v1=(-1,1,0,-1), v2=(-1,0,1,-2)
A=[[1,0],[2,1],[3,2],[1,1]]
v1=[-1,1,0,-1]; v2=[-1,0,1,-2]
print("check v1.A, v2.A:", sp.Matrix([v1])*sp.Matrix(A), sp.Matrix([v2])*sp.Matrix(A))
K=W.int_left_kernel(A)
print("kernel basis:", K)
print("snf invariants (saturation cert):", W.snf_invariants(K))
# non-saturated input to saturate_rows: rows 2*v1, v1+2*v2-ish
B=[[-2,2,0,-2],[-3,1,2,-5]]  # 2*v1, v1+2*v2 (non-saturated)
Bs=W.saturate_rows(B)
print("saturated:", Bs, "snf:", W.snf_invariants(Bs))
qs=sp.Matrix(Bs).nullspace()
print("sat rows still in kernel of A:", sp.Matrix(Bs)*sp.Matrix(A) == sp.zeros(len(Bs),2))

print("=== planted LLL relation at 500 digits (gram='exact') ===")
mp.mp.dps=540
xs=[mp.pi, mp.e, mp.log(2), mp.zeta(3), mp.sqrt(2)]
rel=[2,-3,5,7,-11]
y=-sum(r*x for r,x in zip(rel,xs))/13
vals=xs+[y]
n=len(vals); P=460
C=mp.power(10,P)
B=[]
for i,v in enumerate(vals):
    row=[0]*(n+1); row[i]=1; row[n]=int(mp.nint(v*C)); B.append(row)
import flint
M=flint.fmpz_mat(B)
t0=time.time(); R=M.lll(gram='exact'); t1=time.time()
row0=[int(R[0,j]) for j in range(n)]
print("first LLL row coeffs:", row0, "expected +-", rel+[13], f"({t1-t0:.2f}s)")
res=sum(row0[i]*vals[i] for i in range(n))
print("residual:", mp.nstr(res,5))
ok = all(row0[i]==-s for i,s in zip(rel+[13], row0)) or all(row0[i]==s for i,s in zip(rel+[13],row0))
print("planted relation recovered:", ok)
# pslq cross-check
rq=mp.pslq(vals, tol=mp.mpf(10)**(-P+30), maxcoeff=10**6)
print("pslq:", rq)

print("=== L(chi8,5): Hurwitz vs nsum vs M4i value ===")
P=300
v1=W.Lchi_hurwitz(8,P); v2=W.Lchi_nsum(8,P)
mp.mp.dps=60
print("hurwitz:", mp.nstr(v1,40))
print("diff hurwitz-nsum:", mp.nstr(v1-v2,5))
import json
ref=json.load(open("../M4i_w5_noL/Lchi8_5.json"))["value"]
mp.mp.dps=80
print("match M4i first 60 digits:", mp.nstr(v1,62)[:62]==ref[:62])
print("M4i head:", ref[:70])
print("mine head:", mp.nstr(v1,72))
