import mpmath, pickle
from mpmath import mp, mpf, log, exp, quad, nstr
mp.dps = 130
g1 = lambda s: log(1+exp(-2*s))
g1p = lambda s: -2*exp(-2*s)/(1+exp(-2*s))
def Iab(a,b): return quad(lambda s: s**a*g1(s)**b, [0,mpf('0.5'),1,2,4,8,16,40])
def K2int(m): return quad(lambda s: s**m*g1(s)*g1p(s)**2, [0,mpf('0.5'),1,2,4,8,16,40])
def K3int(m): return quad(lambda s: s**m*g1(s)**2*g1p(s)**2, [0,mpf('0.5'),1,2,4,8,16,40])
needI = [(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(0,8),(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),
         (2,2),(2,3),(2,4),(2,5),(2,6),(3,2),(3,3),(3,4),(3,5),(4,2),(4,3),(4,4),(5,2),(5,3),(6,2)]
I = {}
for key in needI:
    I[key] = Iab(*key)
K2 = {m: K2int(m) for m in [6,7]}
K3 = {m: K3int(m) for m in [6]}
# self-check: I(0,2) = zeta(3)/8
from mpmath import zeta
print("I(0,2)-zeta3/8 =", nstr(I[(0,2)]-zeta(3)/8,5))
pickle.dump({'I':{str(k):nstr(v,120) for k,v in I.items()},
             'K2':{str(k):nstr(v,120) for k,v in K2.items()},
             'K3':{str(k):nstr(v,120) for k,v in K3.items()}}, open('ckpt_integrals130.pkl','wb'))
print("saved", len(I), "integrals")
