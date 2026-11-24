import mpmath as mp, json
from mpmath import mpf
from fractions import Fraction as Fr
mp.mp.dps=75
WORK='research/08_theory/K4d_identify'
import sys; sys.path.insert(0,WORK)
from pslq_lib import pslq
z3=mp.zeta(3); pi=mp.pi; l2=mp.log(2); l3=mp.log(3); sq5=mp.sqrt(5); l43=mp.log(mpf(4)/3)
def sermul(a,b,K):
    o=[mpf(0)]*(K+1)
    for i in range(K+1):
        for j in range(i+1): o[i]+=a[j]*b[i-j]
    return o
def serinv(a,K):
    o=[mpf(0)]*(K+1); o[0]=1/a[0]
    for n in range(1,K+1): o[n]=-sum(a[j]*o[n-j] for j in range(1,n+1))/a[0]
    return o
L=(mpf(-6)/5)*l2+(mpf(13)/5)*l3
u=[mpf(0),1,-1,1,-1];Lu=[L*x for x in u];e=[mpf(0)]*5;e[0]=1;pw=[mpf(0)]*5;pw[0]=1
for k in range(1,5):
    pw=sermul(pw,Lu,4)
    for i in range(4,-1,-1): e[i]+=pw[i]/mp.factorial(k)
finv=serinv(e,4)
def getT(tag):
    SA=[mpf(s) for s in json.load(open(WORK+'/modelA_%s.json'%tag))['Sdk']]
    return sermul(list(SA),finv,4)
T192=getT('h192')
t2=T192[2]/sq5
Li2m13=mp.polylog(2,mpf(-1)/3)
# VERIFY T2 closed form
t2_form = Fr(-254,675) + Fr(79,8100)*pi**2 + Fr(-118,75)*Li2m13 + Fr(-4,15)*l43 + Fr(76,75)*l43**2
t2_form = mpf(str(float(t2_form))) if isinstance(t2_form,Fr) else t2_form
# build properly with mpf
t2_form = mpf(Fr(-254,675)) + mpf(Fr(79,8100))*pi**2 + mpf(Fr(-118,75))*Li2m13 + mpf(Fr(-4,15))*l43 + mpf(Fr(76,75))*l43**2
print("T2/√5 h192   =", mp.nstr(t2,58))
print("T2/√5 closed =", mp.nstr(t2_form,58))
print("  |diff| =", mp.nstr(abs(t2-t2_form),3))
# So T2 = sqrt5 * that; and T2/S0 = T2/(2sqrt5/3) = (3/(2 sqrt5)) T2 = (3/2)*t2_form... wait
# T2/S0 = T2/(2sqrt5/3) = (3 T2)/(2 sqrt5) = (3/2)*(T2/sqrt5) = (3/2)*t2_form
print("\nT2/S0 = (3/2)*t2/√5:")
print("  T2/S0 closed =", mp.nstr(mpf(3)/2*t2_form,40))
