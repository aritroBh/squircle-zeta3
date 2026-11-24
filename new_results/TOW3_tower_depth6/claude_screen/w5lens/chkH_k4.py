import pickle, sys
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps=30
R='<project>/research/08_theory/'
sys.path.insert(0,R+'WTS_weight_k_structure/theory')
K={k:pickle.load(open(R+f'WTS_weight_k_structure/theory/kern_{k}.pkl','rb')) for k in range(4,10)}
# kernel-level: is K4 on the cubic through 5..8 ?
POLES=[(Fr(0),Fr(0)),(Fr(1),Fr(1)),(Fr(-1),Fr(-1)),(Fr(1),Fr(-1))]
off=0
for key in [(4,0),(3,1),(2,2)]:
    for pole in POLES:
        for c in range(2):
            g=lambda k:K[k][key]['real'].get(pole,(Fr(0),)*4)[c]
            pred=4*g(5)-6*g(6)+4*g(7)-g(8)   # cubic through 5..8 at k=4
            if pred!=g(4): off+=1
print('kernel components off the cubic at k=4:',off)
s2=mp.sqrt(2); a=1+s2; b=s2-1
t=lambda v:(v*v+2*v-1)/(4*v)
def W5(k):
    tot=mp.mpf(0)
    for key in [(4,0),(3,1),(2,2)]:
        i,j=key; P=K[k][key]['real']
        def f(v):
            tv=t(v); kern=mp.mpf(0)
            for (r0,r1),co in P.items():
                r=mp.mpf(r0.numerator)/r0.denominator+(mp.mpf(r1.numerator)/r1.denominator)*s2
                c=(mp.mpf(co[0].numerator)/co[0].denominator)+(mp.mpf(co[1].numerator)/co[1].denominator)*s2
                kern+=c/(v-r)
            return kern*mp.log(tv)**i*(mp.log(1-tv)**j if j else 1)
        tot+=mp.quad(f,[b,1,a])
    return tot.real
W={k:W5(k) for k in range(4,9)}
P4=4*W[5]-6*W[6]+4*W[7]-W[8]
print('W5(4) from kern_4 =',mp.nstr(W[4],20),' cubic extrapolation P(4) =',mp.nstr(P4,20),' P(4)-W5(4) =',mp.nstr(P4-W[4],20))
