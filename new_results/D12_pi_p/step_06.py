# step_06.py -- TASK 2/4: PSLQ identifications (95 digits) and sqrt2-absence test
import mpmath as mp, pickle
from math import gcd
mp.mp.dps = 95
N = pickle.load(open('/tmp/Nvals.pkl','rb'))['N']
pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3); z5=mp.zeta(5)
Li4=mp.polylog(4,mp.mpf('0.5')); Li5=mp.polylog(5,mp.mpf('0.5'))
def identify(tag,t,b,nm):
    rel=mp.pslq([t]+b,maxcoeff=10**10,maxsteps=10000)
    g=0
    for x in rel: g=gcd(g,x)
    rel=[x//g for x in rel]
    val=-sum(rel[i+1]*b[i] for i in range(len(b)))/rel[0]
    print(tag, rel, "residual=", mp.nstr(abs(val-t),2))
identify('N2',N[2],[L2**3,pi**2],[])
identify('N3',N[3],[L2**4,pi**2,pi**2*L2,z3],[])
identify('N4',N[4],[L2**5,pi**2,pi**2*L2,pi**2*L2**2,z3,z3*L2,pi**4],[])
# J-evaluations at weights 4 and 5 (Euler sums)
def J(m,k): return mp.quad(lambda s: s**k*mp.log(1+mp.e**(-s))**m,[0,mp.inf])
W4=[L2**4,pi**2*L2**2,pi**4,z3*L2,Li4]
identify('J30',J(3,0),W4,[]); identify('J21',J(2,1),W4,[])
W5=[L2**5,pi**2*L2**3,pi**4*L2,z3*L2**2,z3*pi**2,z5,Li4*L2,Li5]
identify('J40',J(4,0),W5,[]); identify('J31',J(3,1),W5,[]); identify('J22',J(2,2),W5,[])
# N5 closed form check
N5form=(L2**6/720 + pi**2/24*(1-L2+L2**2/2-L2**3/6) + z3/4*(1-L2+L2**2/2)
        + pi**4/2880*(1-L2) - 7*pi**4/1920 - z3*pi**2/48 + z5/8)
print("N5 form residual:", mp.nstr(abs(N5form-N[5]),2))
# sqrt2 absence: give PSLQ BOTH alphabets
sq2=mp.sqrt(2); las=mp.log(1+sq2); li2m=mp.polylog(2,sq2-1)
bE=[sq2,las,sq2*las,li2m,las**2,L2**4,pi**2,pi**2*L2,z3]
rel=mp.pslq([N[3]]+bE,maxcoeff=10**9,maxsteps=10000)
print("N3 relation in joint basis (first 5 slots = Euclidean constants):", rel)
