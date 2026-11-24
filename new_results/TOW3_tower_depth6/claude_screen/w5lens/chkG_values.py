# chkG_values.py -- numerical cross-check: W5(k) from exact kernels by quadrature (k=5..9), finite differences E1,E2,E3, vs paper's constants
import pickle, sys
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps=40
R='<project>/research/08_theory/'
sys.path.insert(0,R+'WTS_weight_k_structure/theory')
s2=mp.sqrt(2); a=1+s2; b=s2-1
t=lambda v:(v*v+2*v-1)/(4*v)
KS={k:pickle.load(open(R+f'WTS_weight_k_structure/theory/kern_{k}.pkl','rb')) for k in range(5,10)}
def W5(k):
    tot=mp.mpf(0)
    for key in [(4,0),(3,1),(2,2)]:
        i,j=key
        P=KS[k][key]['real']
        def f(v):
            tv=t(v)
            kern=mp.mpf(0)
            for (r0,r1),co in P.items():
                r=mp.mpf(r0.numerator)/r0.denominator+(mp.mpf(r1.numerator)/r1.denominator)*s2
                c=(mp.mpf(co[0].numerator)/co[0].denominator)+(mp.mpf(co[1].numerator)/co[1].denominator)*s2
                kern+=c/(v-r)
            return kern*mp.log(tv)**i*(mp.log(1-tv)**j if j else 1)
        tot+=mp.quad(f,[b,1,a])
    return tot
W={k:W5(k) for k in range(5,10)}
for k in W: print('W5(%d) = %s'%(k,mp.nstr(W[k],30)))
d1=W[6]-W[5]; d2=W[7]-2*W[6]+W[5]; d3=W[8]-3*W[7]+3*W[6]-W[5]
d4=W[9]-4*W[8]+6*W[7]-4*W[6]+W[5]
print('E1 =',mp.nstr(d1,25)); print('E2 =',mp.nstr(d2,25)); print('E3 =',mp.nstr(d3,25)); print('4th diff =',mp.nstr(d4,5))
lam=mp.log(2)
E3c=lam**5/48+5*mp.pi**2*lam**3/144+mp.pi**4*lam/64+5*mp.zeta(3)*lam**2/16+5*mp.pi**2*mp.zeta(3)/96+15*mp.zeta(5)/32
print('paper E3 closed form =',mp.nstr(E3c,25),' diff =',mp.nstr(d3-E3c,5))
print('paper E1 = 8.36025639176937361007  diff', mp.nstr(d1-mp.mpf('8.36025639176937361007'),5))
print('paper E2 = 8.70728401062393015940  diff', mp.nstr(d2-mp.mpf('8.70728401062393015940'),5))
