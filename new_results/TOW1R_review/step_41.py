import subprocess, os
work = "<project>/research/08_theory/TOW1R_review/"
step = r'''
import sys, os
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from myengine import *
import sympy as sp
i=sp.I
t = { (1,0): (-1-i)/2, (2,0): (-1+i)/2, (1,1): sp.Rational(3,4)*(1-i),
      (2,1): -sp.Rational(3,4)-sp.Rational(3,4)*i, (1,2): 1+sp.Rational(3,4)*i, (2,2): 1-sp.Rational(3,4)*i }
def tc1c(n): return (-1)**n*t[(2,n)]
def rfun(dd,m,conjg=False):
    s=-2*dd if not conjg else 2*dd
    return 2*(dd-2*m)+s*i-4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd,m,conjg=False):
    s=-1 if not conjg else 1
    return sp.Rational(1,2)*((3*dd-3)/2+s*i*rfun(dd,m,conjg)/2-(dd-2)**2)
def proveE3_feed(P, delta):
    p=sp.Symbol('p'); rB0=2*i*(4*p-1)/(p+1)
    def erho(dl,conjg=False):
        r0=rB0 if not conjg else sp.conjugate(rB0)
        if dl==0: return r0
        if dl==1: return (-1 if not conjg else 1)*i*(2*p-1)*r0
        if dl==2: return R2(2*p+1,p+1,conjg)*r0
    def rho_pp(dl,conjg=False):
        if dl==0: return sp.Integer(1)
        if dl==1: return (-1 if not conjg else 1)*i*(2*p-2)
        if dl==2: return R2(2*p,p,conjg)
    def ehat(dl,s,conjg=False):
        Pp=p+1
        if (dl,s)==(0,0): return erho(0,conjg)*(t[(2,0)] if conjg else t[(1,0)])*Pp/(4*p-1)
        if (dl,s)==(1,0): return (erho(1,conjg)*(t[(2,0)] if conjg else t[(1,0)])*Pp/(4*p-2)+erho(0,conjg)*((-t[(2,1)]) if conjg else t[(1,1)])*Pp/(4*p-1))
        if (dl,s)==(2,0): return (erho(2,conjg)*(t[(2,0)] if conjg else t[(1,0)])*Pp/(4*p-3)+erho(1,conjg)*((-t[(2,1)]) if conjg else t[(1,1)])*Pp/(4*p-2)+erho(0,conjg)*(t[(2,2)] if conjg else t[(1,2)])*Pp/(4*p-1))
        if (dl,s)==(2,1): return (-erho(0,conjg)/(4*p-1))*(t[(2,0)] if conjg else t[(1,0)])*Pp/(4*p-3)
    def Delta(dl,conjg=False):
        return sp.expand(rho_pp(dl,conjg)-sum((-1)**s*ehat(dl,s,conjg) for s in range(0,dl//2+1)))
    def Dc(dl): return -Delta(dl,conjg=True)
    def gfun(mu):
        terms={4*p-4:[(0,0)],4*p-5:[(0,1),(1,0)],4*p-6:[(0,2),(1,1),(2,0)]}
        return sum(-Dc(dp)*tc1c(n)/(4*p-3-dp) for (dp,n) in terms[mu])
    def feed(dl):
        direct=p*sum(t[(2,dl-dp)]*Delta(dp)/(4*p-3-dp) for dp in range(0,dl+1))
        if dl==0: ft=-p*(gfun(4*p-4)*(-1))
        elif dl==1: ft=-p*(gfun(4*p-5)*(-1)+gfun(4*p-4)*i*(4*p-6))
        elif dl==2: ft=-p*(gfun(4*p-6)*(-1)+gfun(4*p-5)*i*(4*p-7)+gfun(4*p-4)*sp.binomial(4*p-6,2))
        return sp.expand(direct+ft)
    return sp.simplify(feed(delta).subs(p,P).replace(sp.binomial(4*P-6,2),(4*P-6)*(4*P-7)/2))
# my numeric Delta and feed
from math import comb as C
def DeltaN(p,delta):
    out=rho(p,p,delta)
    for s in range(0,delta//2+1): out=out-ehat(p,delta,s).scale((-1)**s)
    return out
def FoldN(p,delta):
    m=4*p-4-delta; tot=ZERO
    for mu in range(m,4*p-3):
        g=ZERO
        for dp in range(0,4*p-4-mu+1):
            g=g-DeltaN(p,dp)*tau1(4*p-4-mu-dp).scale(Fr(1,4*p-3-dp))
        tot=tot+g.conj_i()*(IQ**(-mu))*Fr(C(mu-2,mu-m))*(IQ**(m-2))
    return tot.scale(-p)
def feedN(p,delta):
    direct=ZERO
    for dp in range(0,delta+1): direct=direct+DeltaN(p,dp)*tau2(delta-dp).scale(Fr(p,4*p-3-dp))
    return direct+FoldN(p,delta)
for P in (3,4,5):
    for delta in (0,1,2):
        pv = proveE3_feed(P,delta)
        mine = (feedN(P,delta)/rho(P,P,0)).sympy()
        print(f"E3 p={P} delta={delta}: prove-feed == mine/rho : {sp.simplify(sp.expand(pv)-sp.expand(mine))==0}  (value={sp.simplify(pv)})")
'''
with open(work+"step_14.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_14.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=600)
print(r.stdout); print("STDERR:",r.stderr[:2500])
with open(work+"progress.log","a") as fh: fh.write("step_14: cross-check prove_E3 feed vs my numeric (normalized), p=3,4,5\n")
