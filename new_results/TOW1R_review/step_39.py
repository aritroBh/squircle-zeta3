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
def tc1c(n): return (-1)**n * t[(2,n)]
def rfun(dd,m,conjg=False):
    s=-2*dd if not conjg else 2*dd
    return 2*(dd-2*m)+s*i-4*(dd-1)*(dd-2*m)/(2*dd-3)
def R2(dd,m,conjg=False):
    s=-1 if not conjg else 1
    return sp.Rational(1,2)*((3*dd-3)/2 + s*i*rfun(dd,m,conjg)/2 - (dd-2)**2)

# ---- E2b non-adjacent ----
def E2bnonadj_LR(DD,PP):
    d,p=sp.symbols('d p')
    LHS=R2(d,p)-1/(2*d-3)
    rA0=i*(2*d-1)/(p+1); rB0=-i*(2*d-1)/(d-p+1)
    RHS=((p+1)*(t[(1,2)]*rA0/(2*d-1)+t[(1,1)]*(-i*(d-1))*rA0/(2*d-2)+t[(1,0)]*(R2(d+1,p+1)*rA0)/(2*d-3))
        +(d-p+1)*(t[(2,2)]*rB0/(2*d-1)+t[(2,1)]*(-i*(d-1))*rB0/(2*d-2)+t[(2,0)]*(R2(d+1,p)*rB0)/(2*d-3)))
    return sp.simplify(LHS.subs({d:DD,p:PP})), sp.simplify(RHS.subs({d:DD,p:PP}))
for (DD,PP) in [(5,4),(6,5),(6,4),(7,6)]:
    QQ=DD-PP
    L,R=E2bnonadj_LR(DD,PP)
    mineL=((e_off(PP,QQ,2,0)+e_off(PP,QQ,2,1))/rho(PP,QQ,0)).sympy()
    mineR=(LambdaE(PP,QQ,2,0)/rho(PP,QQ,0)).sympy()
    print(f"E2b-nonadj d={DD} p={PP}: LHS==mine {sp.simplify(L-mineL)==0}, RHS==mine {sp.simplify(sp.expand(R)-sp.expand(mineR))==0}")

# ---- E2b adjacent ----
def E2badj_LR(P):
    p=sp.Symbol('p'); d=2*p-1
    rat={(p+1,p-1):i*(2*d-1)/(p+1),(p+1,p):(2*d-1)*(2*d+1)/((p+1)*p)}
    def e_ratio(key,delta,conjg=False):
        D=key[0]+key[1]; r0=rat[key] if not conjg else sp.conjugate(rat[key])
        if delta==0: return r0
        if delta==1: return (-1 if not conjg else 1)*i*(D-2)*r0
        if delta==2: return R2(D,key[0],conjg)*r0
    def e_off_slope_ratio(conjg=False):
        return -rat[(p+1,p)]/(2*(d+2)-3) if not conjg else -sp.conjugate(rat[(p+1,p)])/(2*(d+2)-3)
    def ehat(delta,s,conjg=False):
        Pp=p+1
        if (delta,s)==(0,0): return e_ratio((Pp,p),0,conjg)*(t[(2,0)] if conjg else t[(1,0)])*Pp/(4*p-1)
        if (delta,s)==(1,0): return (e_ratio((Pp,p),1,conjg)*(t[(2,0)] if conjg else t[(1,0)])*Pp/(4*p-2)+e_ratio((Pp,p),0,conjg)*((-t[(2,1)]) if conjg else t[(1,1)])*Pp/(4*p-1))
        if (delta,s)==(2,0): return (e_ratio((Pp,p),2,conjg)*(t[(2,0)] if conjg else t[(1,0)])*Pp/(4*p-3)+e_ratio((Pp,p),1,conjg)*((-t[(2,1)]) if conjg else t[(1,1)])*Pp/(4*p-2)+e_ratio((Pp,p),0,conjg)*(t[(2,2)] if conjg else t[(1,2)])*Pp/(4*p-1))
        if (delta,s)==(2,1): return e_off_slope_ratio(conjg)*(t[(2,0)] if conjg else t[(1,0)])*Pp/(4*p-3)
    def ed(delta,conjg=False):
        if delta in (0,1): return ehat(delta,0,conjg)
        if delta==2: return ehat(2,0,conjg)-ehat(2,1,conjg)
    LHS=R2(d,p)-1/(2*d-3)
    feed_DL1=(p+1)*(t[(1,2)]*rat[(p+1,p-1)]/(2*d-1)+t[(1,1)]*e_ratio((p+1,p-1),1)/(2*d-2)+t[(1,0)]*e_ratio((p+1,p-1),2)/(2*d-3))
    feed_DL2=p*(t[(2,2)]*ed(0)/(2*d-1)+t[(2,1)]*ed(1)/(2*d-2)+t[(2,0)]*ed(2)/(2*d-3))
    g_m=-(ed(0,True)*t[(2,2)]/(4*p-3)+ed(1,True)*(-t[(2,1)])/(4*p-4)+ed(2,True)*t[(2,0)]/(4*p-5))
    g_m1=-(ed(0,True)*(-t[(2,1)])/(4*p-3)+ed(1,True)*t[(2,0)]/(4*p-4))
    g_m2=-(ed(0,True)*t[(2,0)]/(4*p-3))
    feed_fold=p*(g_m-i*(4*p-7)*g_m1-sp.binomial(4*p-6,2)*g_m2)
    RHS=feed_DL1+feed_DL2+feed_fold
    return sp.simplify(LHS.subs(p,P)), sp.simplify(RHS.subs(p,P).replace(sp.binomial(4*P-6,2),(4*P-6)*(4*P-7)/2))
for P in (3,4,5):
    L,R=E2badj_LR(P)
    mineL=((e_off(P,P-1,2,0)+e_off(P,P-1,2,1))/rho(P,P-1,0)).sympy()
    mineR=(LambdaE(P,P-1,2,0)/rho(P,P-1,0)).sympy()
    print(f"E2b-adj p={P}: LHS==mine {sp.simplify(L-mineL)==0}, RHS==mine {sp.simplify(sp.expand(R)-sp.expand(mineR))==0}")
'''
with open(work+"step_12.py","w") as f: f.write(step)
env=dict(os.environ); env["PYTHONDONTWRITEBYTECODE"]="1"
r=subprocess.run(["nice","-n","19","<project>/.venv/bin/python","step_12.py"],
                 cwd=work,capture_output=True,text=True,env=env,timeout=900)
print(r.stdout); print("STDERR:",r.stderr[:2500])
with open(work+"progress.log","a") as fh: fh.write("step_12: cross-check E2b non-adj + adj scripts vs my numeric LambdaE\n")
