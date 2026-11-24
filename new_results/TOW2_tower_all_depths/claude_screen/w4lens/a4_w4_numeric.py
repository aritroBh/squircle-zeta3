# Item (4) cross-check: assemble W4(k)=sum int K^{(k)}_{ij} l1^i l2^j dv from exact kernel data (mpmath, 40 digits),
# compare with the paper's W4''+(k-4)D1+C(k-4,2)D2 (closed forms D1,D2 evaluated independently), k=3..14.
import pickle, os, sys
from fractions import Fraction as Fr
from mpmath import mp, mpf, sqrt, log, quad, polylog, zeta, pi, binomial
sys.dont_write_bytecode = True
mp.dps = 40
TH="<project>/research/08_theory/WTS_weight_k_structure/theory"
RD="<project>/research/08_theory/TWOEND1_weight_poles/red"
def path(k):
    for D in (TH,RD):
        p=os.path.join(D,f"kern_{k}.pkl")
        if os.path.exists(p): return p
r2=sqrt(2); a=1+r2; b=r2-1
def l1(v): return log(v-b)+log(v+a)-log(4*v)
def l2(v): return log(a-v)+log(v+b)-log(4*v)
def pole_val(p): return mpf(p[0].numerator)/p[0].denominator + (mpf(p[1].numerator)/p[1].denominator)*r2
I={}
def Iint(i,j,p):
    key=(i,j,p)
    if key not in I:
        r=pole_val(p)
        f=lambda v: l1(v)**i*l2(v)**j/(v-r)
        I[key]=quad(f,[b,1])+quad(f,[1,a])
    return I[key]
def W4(k):
    K=pickle.load(open(path(k),'rb'))
    tot=mpf(0)
    for key in [(3,0),(2,1)]:
        for p,c in K[key]['real'].items():
            assert c[2]==0 and c[3]==0
            if c[0]==0 and c[1]==0: continue
            cv=mpf(c[0].numerator)/c[0].denominator+(mpf(c[1].numerator)/c[1].denominator)*r2
            tot+=cv*Iint(key[0],key[1],p)
    return tot
lam=log(2); L=log(a)
z3=zeta(3); b_=b
D2=-lam**4/24-pi**2*lam**2/24-pi**4/160-z3*lam/4
D1=(-mpf(5)/48*lam**4-lam**3*L/12-lam**2*L**2/4+lam*L**3/6+L**4/3-mpf(7)/96*pi**2*lam**2
    +mpf(5)/24*pi**2*lam*L-mpf(5)/24*pi**2*L**2-pi**4/192-mpf(3)/8*z3*lam-mpf(7)/16*z3*L
    -lam*polylog(3,b_)-2*L*polylog(3,b_)-lam*polylog(3,-b_)+2*L*polylog(3,-b_)+4*L*polylog(3,1/r2))
Wd=W4(4)
print("W4(4) from kernel data =",Wd)
print("paper W4''            = 1.42445367600887292151442249125727")
res={}
for k in range(3,13):
    res[k]=W4(k)
    form=Wd+(k-4)*D1+binomial(k-4,2)*D2
    print(f"k={k:2d}  W4(k)={mp.nstr(res[k],30)}  formula={mp.nstr(form,30)}  diff={mp.nstr(res[k]-form,5)}")
print("D1(num) =",mp.nstr(D1,30)," from data W4(5)-W4(4) =",mp.nstr(res[5]-res[4],30))
print("D2(num) =",mp.nstr(D2,30)," from data Delta^2     =",mp.nstr(res[6]-2*res[5]+res[4],30))
print("offset at k=3: W4''-D1+D2-W4(3) =",mp.nstr(Wd-D1+D2-res[3],32)," (paper 0.15446753106687400637871015148911)")
