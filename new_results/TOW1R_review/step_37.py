import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode=True
WTS="<project>/research/08_theory/WTS_weight_k_structure/theory/"
sys.path.insert(0,WTS); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import comb, factorial
from fastred import *
import sympy as sp
i=sp.I
t = { (1,0): (-1-i)/2, (2,0): (-1+i)/2, (1,1): sp.Rational(3,4)*(1-i),
      (2,1): -sp.Rational(3,4)-sp.Rational(3,4)*i, (1,2): 1+sp.Rational(3,4)*i, (2,2): 1-sp.Rational(3,4)*i }
# --- my numeric engine ---
def doublefact(n):
    r=1
    for x in range(n,0,-2): r*=x
    return r
def Kd(dd): return (QSI(Fr((-1)**dd*doublefact(2*dd-3),2*factorial(dd)))*(IQ**(dd+1)))*QSI(0,Fr(1,2))
def rp(dd,pp): return QSI(2*(dd-2*pp),0,-2*dd,0)-QSI(Fr(4*(dd-1)*(dd-2*pp),2*dd-3))
def rho(pp,qq,delta):
    dd=pp+qq; mult=2 if pp>qq else 1
    r0=Kd(dd)*QSI(Fr((-1)**pp*comb(dd,pp)*mult))
    if delta==0: return r0
    if delta==1: return r0*QSI(0,0,-(dd-2),0)
    val=QSI(Fr(3*dd-3,2))-rp(dd,pp)*(IQ.scale(Fr(1,2)))-QSI((dd-2)**2)
    return (r0*val).scale(Fr(1,2))
def tau1(n):
    s=n+1; return ((-1)**s)*((IQ-BQ)**(-(s+1))+(IQ+AQ)**(-(s+1))-IQ**(-(s+1)))
def tau2(n):
    s=n+1; return ((-1)**s)*((IQ-AQ)**(-(s+1))+(IQ+BQ)**(-(s+1))-IQ**(-(s+1)))
def e_off(pp,qq,delta,s):
    if (delta,s)==(0,0): return rho(pp,qq,0)
    if (delta,s)==(1,0): return rho(pp,qq,1)
    if (delta,s)==(2,0): return rho(pp,qq,2)
    if (delta,s)==(2,1): return rho(pp,qq,0).scale(Fr(-1,2*(pp+qq)-3))
    return ZERO
def ehat(pp,delta,s):
    P=pp+1
    if (delta,s)==(0,0): return e_off(pp+1,pp,0,0)*tau1(0).scale(Fr(P,4*pp-1))
    if (delta,s)==(1,0): return (e_off(pp+1,pp,1,0)*tau1(0).scale(Fr(P,4*pp-2))+e_off(pp+1,pp,0,0)*tau1(1).scale(Fr(P,4*pp-1)))
    if (delta,s)==(2,0): return (e_off(pp+1,pp,2,0)*tau1(0).scale(Fr(P,4*pp-3))+e_off(pp+1,pp,1,0)*tau1(1).scale(Fr(P,4*pp-2))+e_off(pp+1,pp,0,0)*tau1(2).scale(Fr(P,4*pp-1)))
    if (delta,s)==(2,1): return e_off(pp+1,pp,2,1)*tau1(0).scale(Fr(P,4*pp-3))
    return ZERO
def e_any(pp,qq,delta,s): return ehat(pp,delta,s) if pp==qq else e_off(pp,qq,delta,s)
def Ediag_conv(pp,qq,delta,sigma):
    if pp!=qq: return e_off(pp,qq,delta,sigma)
    out=ZERO; tt=0
    while sigma+tt<=delta//2:
        out=out+ehat(pp,delta,sigma+tt).scale((-1)**tt); tt+=1
    return out
def LambdaE(pp,qq,delta,sigma):
    dd=pp+qq; tot=ZERO
    for dp in range(0,delta+1):
        tot=tot+e_any(pp+1,qq,dp,sigma)*tau1(delta-dp).scale(Fr(pp+1,2*dd-dp-1))
        if pp>qq:
            tot=tot+Ediag_conv(pp,qq+1,dp,sigma)*tau2(delta-dp).scale(Fr(qq+1,2*dd-dp-1))
    if qq==pp-1:
        m=4*pp-4-delta; ft=ZERO
        for mu in range(m,4*pp-3):
            g=ZERO
            for dp in range(0,4*pp-4-mu+1):
                g=g-Ediag_conv(pp,pp,dp,sigma)*tau1(4*pp-4-mu-dp).scale(Fr(1,4*pp-3-dp))
            ft=ft+g.conj_i()*(IQ**(-mu))*Fr(comb(mu-2,mu-m))*(IQ**(m-2))
        tot=tot+ft.scale(-pp)
    return tot

# --- prove_E1_adjacent symbolic LHS,RHS ---
def proveE1_LR(P):
    p=sp.Symbol('p'); dd=2*p-1
    LHS=-i*(dd-2)
    rA0=i*(2*dd-1)/(p+1); rA1=-i*(dd-1)*rA0
    rB0=(2*dd-1)*(2*dd+1)/((p+1)*p); rB1=-i*dd*rB0
    ehat00=t[(1,0)]*(4*p-3)/p
    ehat10=(p+1)*(t[(1,0)]*rB1/(4*p-2)+t[(1,1)]*rB0/(4*p-1))
    feed_DL1=(p+1)*(t[(1,0)]*rA1/(2*dd-2)+t[(1,1)]*rA0/(2*dd-1))
    feed_DL2=p*(t[(2,0)]*ehat10/(2*dd-2)+t[(2,1)]*ehat00/(2*dd-1))
    ehat00c=sp.conjugate(ehat00)
    ehat10c=(p+1)*(t[(2,0)]*sp.conjugate(rB1)/(4*p-2)+(-t[(2,1)])*sp.conjugate(rB0)/(4*p-1))
    g_m1c=-(ehat00c*(-t[(2,1)])/(4*p-3)+ehat10c*t[(2,0)]/(4*p-4))
    g_m0c=-(ehat00c*t[(2,0)]/(4*p-3))
    feed_fold=p*(g_m1c-i*(4*p-6)*g_m0c)
    RHS=feed_DL1+feed_DL2+feed_fold
    return sp.simplify(LHS.subs(p,P)), sp.simplify(RHS.subs(p,P))
for P in (3,4,5):
    L,R = proveE1_LR(P)
    mineL = (e_off(P,P-1,1,0)/rho(P,P-1,0)).sympy()
    mineR = (LambdaE(P,P-1,1,0)/rho(P,P-1,0)).sympy()
    print(f"E1-adj p={P}: proveLHS==mine {sp.simplify(L-mineL)==0}, proveRHS==mine {sp.simplify(sp.expand(R)-sp.expand(mineR))==0}")
'''
with open(work+"step_11.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_11.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=600)
print(r.stdout); print("STDERR:",r.stderr[:2000])
with open(work+"progress.log","a") as fh: fh.write("step_11: cross-check prove_E1_adjacent LHS/RHS == my numeric LambdaE (normalized), p=3,4,5\n")
