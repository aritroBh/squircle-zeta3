import sys, subprocess, re
lib = open('w5lib.py').read()
# replace the broken cobracket function with a correct one
start = lib.index("# ---------------- cobracket (3,2) ----------------")
new_tail = r'''# ---------------- cobracket (3,2) ----------------
def _perm_sign(seq):
    """sign of the permutation that sorts seq (seq must have distinct entries)."""
    inv = 0
    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            if seq[i] > seq[j]: inv += 1
    return -1 if inv % 2 else 1

def cobracket32(T, d, mode="full"):
    """T: dict (i0..i4) -> coeff (letters = formal basis of V, dim d).
    Returns dict ((triple),(pair)) -> coeff for delta_{3,2}: V^ot5 -> Lambda^3 V (x) Lambda^2 V.
    mode 'full': sum over all C(5,3)=10 (3+2)-shuffles with signs (Lie cobracket);
    mode 'cut' : plain deconcatenation (0,1,2)|(3,4) followed by antisymmetrisation."""
    out = {}
    for tup, c in T.items():
        if c == 0: continue
        if mode == "cut":
            splits = [((0,1,2),(3,4),1)]
        else:
            splits = []
            for I in itertools.combinations(range(5),3):
                J = tuple(j for j in range(5) if j not in I)
                sgn = _perm_sign(list(I)+list(J))
                splits.append((I,J,sgn))
        for I,J,sgn in splits:
            tri = [tup[I[0]], tup[I[1]], tup[I[2]]]
            if len(set(tri)) < 3: continue
            pair = [tup[J[0]], tup[J[1]]]
            if pair[0] == pair[1]: continue
            tsgn = _perm_sign(tri)
            psgn = 1 if pair[0] < pair[1] else -1
            key = (tuple(sorted(tri)), tuple(sorted(pair)))
            out[key] = out.get(key, 0) + c*sgn*tsgn*psgn
    return {k:v for k,v in out.items() if v}
'''
open('w5lib.py','w').write(lib[:start] + new_tail)

test = r'''
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import itertools, math, time
import mpmath as mp
import sympy as sp
import w5lib as W
from w5lib import QD

print("=== QD identities ===")
w2=QD(0,1,8); b=w2-1
print("(1+w2)^2 == 3+2w2:", (1+w2)**2 == QD(3,2,8))
phi=QD(Fraction(1,2),Fraction(1,2),5)
print("phi^2 == phi+1:", phi**2 == phi+1)
r3=QD(0,1,12)
print("(1+r3)^2 == 2*(2+r3):", (1+r3)**2 == 2*(2+r3))
print("b^4 == 17-12w2:", b**4 == QD(17,-12,8))
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
        x=QD(1,0,D)
        for g,e in zip(C['letters'],exps): x=x*(g**e)
        if random.choice([0,1]): x=-x
        f=W.factor_sunit(x,C)
        assert f==exps, (D,exps,f)
        # reconstruct |x| from coords
        y=QD(1,0,D)
        for g,e in zip(C['letters'],f): y=y*(g**e)
        assert y==x or y==-x
        ok+=1
    # rejection tests
    if D==8: bad=[QD(3,1,8)]            # N=7
    if D==5: bad=[QD(3,2,5)]            # N=-11
    if D==12: bad=[QD(4,1,12)]          # N=13
    if D==13: bad=[QD(1,1,13)]          # N=-12= -2^2*3 -> 3 not allowed
    for x in bad:
        assert W.factor_sunit(x,C) is None, (D,x)
    print(f"D={D}: 20/20 random S-units factored+reconstructed, rejections ok; {C['Sdesc']}")

print("=== saturation & kernel on planted matrices ===")
# planted: A (4x2) rank2, left kernel spanned by v1=(1,-2,1,0), v2=(0,1,-2,1) but give NON-saturated basis
A=[[1,0],[2,1],[3,2],[1,1]]
# check: (1,-2,1,0).A = (1-4+3, -2+2)=(0,0) ok; (0,1,-2,1).A=(2-6+1,1-4+1)=(-3,-2) NOT kernel. fix:
v2=[ -1,1,1,-3]
print("check v2.A:", sp.Matrix([v2])*sp.Matrix(A))
K=W.int_left_kernel(A)
print("kernel basis:", K)
print("snf invariants (saturation cert):", W.snf_invariants(K))
# non-saturated input to saturate_rows: rows 2*v1, v1+2*v2-ish
B=[[2,-4,2,0],[ -2,2,2,-6]]
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
'''
open('steps/step_06_tests.py','w').write(test)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_06_tests.py'],capture_output=True,text=True,timeout=1200)
print(r.stdout); print("STDERR:", r.stderr[-2000:])
with open('progress.log','a') as f: f.write("step_06_tests: arithmetic/factor/saturation/planted-LLL/L-value tests\n")
