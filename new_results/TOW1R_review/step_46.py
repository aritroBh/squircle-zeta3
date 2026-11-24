import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS)
from fractions import Fraction as Fr
from genintegrand import S_poly
from fastred import *
import sympy as sp
i=sp.I
# Corollary S3: (S_m'/S_m)(t*) = r_m = 2(d-2m) - 2 d i - 4(d-1)(d-2m)/(2d-3)
tstar = sp.Rational(1,2)+i/2
def S_val(d,m,t): 
    P=S_poly(d,m); return sum(Fr(c)*t**s for s,c in enumerate(P))
def Sp_val(d,m,t):
    P=S_poly(d,m); return sum(Fr(s)*Fr(c)*t**(s-1) for s,c in enumerate(P) if s>=1)
# use exact complex arithmetic via QSI at tstar=(1+i)/2
tsQ = QSI(Fr(1,2),0,Fr(1,2),0)
def S_valQ(d,m):
    P=S_poly(d,m); out=ZERO
    for s,c in enumerate(P): out = out + QSI(c)*(tsQ**s)
    return out
def Sp_valQ(d,m):
    P=S_poly(d,m); out=ZERO
    for s,c in enumerate(P):
        if s>=1: out = out + QSI(Fr(s*c))*(tsQ**(s-1))
    return out
def r_m(d,m): return QSI(2*(d-2*m),0,-2*d,0)-QSI(Fr(4*(d-1)*(d-2*m),2*d-3))
bad=0
for d in range(1,15):
    for m in range(0,d+1):
        Sv=S_valQ(d,m); Spv=Sp_valQ(d,m)
        if Sv.is_zero(): continue
        ratio = Spv/Sv
        if ratio != r_m(d,m): bad+=1; print("S3 mismatch",d,m)
print("Corollary S3 (S'/S)(t*)=r_m, d=1..14:", "0 bad" if bad==0 else f"{bad} bad")

# log-derivative identity for A(v): build A symbolically per (d,m), check (log A)'(i) = -i(d-2)
v=sp.symbols('v'); sq2=sp.sqrt(2); a=1+sq2; b=sq2-1
def A_sym(d,m):
    P = S_poly(d,m); D=len(P)-1
    tv = (v**2+2*v-1)/(4*v)
    St = sum(Fr(c)*tv**s for s,c in enumerate(P))*(4*v)**D
    return -(sq2/2)*8**d/4**D * v**(2*d-1-D) * St * (v**2+1)**(2-2*d)/((v**2-a**2)*(v**2-b**2))
bad2=0
for d in range(2,9):
    for m in range(1,d+1):
        A=A_sym(d,m)
        logA1 = sp.simplify(sp.diff(sp.log(A),v).subs(v,i))
        if sp.simplify(logA1 - (-i*(d-2)))!=0:
            bad2+=1; print("(logA)' mismatch",d,m, sp.simplify(logA1))
print("(log A)'(i) = -i(d-2), d=2..8:", "0 bad" if bad2==0 else f"{bad2} bad")
# second derivative: (logA)''(i) + (logA)'(i)^2 = 2 rho(2)/rho(0) = (3d-3)/2 - (i/2) r_m - (d-2)^2
bad3=0
for d in range(2,9):
    for m in range(1,d+1):
        A=A_sym(d,m)
        l1 = sp.diff(sp.log(A),v).subs(v,i); l2 = sp.diff(sp.log(A),v,2).subs(v,i)
        A2overA = sp.simplify(l2 + l1**2)   # A''/A = (logA)'' + (logA)'^2 ; rho(2)=A''(i)/2 => 2 rho2/rho0 = A''/A
        target = sp.simplify((3*d-3)/2 - i/2*( 2*(d-2*m) - 2*d*i - 4*(d-1)*(d-2*m)/(2*d-3) ) - (d-2)**2)
        if sp.simplify(A2overA - target)!=0:
            bad3+=1; print("rho2 mismatch",d,m, sp.simplify(A2overA), target)
print("2*rho(2)/rho(0) = (3d-3)/2-(i/2)r_m-(d-2)^2, d=2..8:", "0 bad" if bad3==0 else f"{bad3} bad")
'''
with open(work+"step_18.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_18.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=900)
print(r.stdout); print("STDERR:",r.stderr[:2500])
with open(work+"progress.log","a") as fh: fh.write("step_18: Corollary S3 + logA derivs (rho1,rho2 derivation) verified symbolically\n")
