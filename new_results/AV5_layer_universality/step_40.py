import mpmath, pickle
from mpmath import mp, mpf, log, exp, quad, zeta, nstr
mp.dps = 135
g1 = lambda s: log(1+exp(-2*s))
g1p = lambda s: -2*exp(-2*s)/(1+exp(-2*s))
PTS = [0, mpf('0.5'), 1, 2, 4, 8, 16, 32, 56, 80]
def Iab(a,b): return quad(lambda s: s**a*g1(s)**b, PTS)
def K2int(m): return quad(lambda s: s**m*g1(s)*g1p(s)**2, PTS)
def K3int(m): return quad(lambda s: s**m*g1(s)**2*g1p(s)**2, PTS)
needI = [(0,2),(0,3),(0,4),(0,5),(0,6),(0,7),(0,8),(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),
         (2,2),(2,3),(2,4),(2,5),(2,6),(3,2),(3,3),(3,4),(3,5),(4,2),(4,3),(4,4),(5,2),(5,3),(6,2)]
I = {key: Iab(*key) for key in needI}
K2 = {m: K2int(m) for m in [6,7]}
K3 = {m: K3int(m) for m in [6]}
print("I(0,2)-zeta3/8 =", nstr(I[(0,2)]-zeta(3)/8,5))
print("I(0,3) =", nstr(I[(0,3)],50))
pickle.dump({'I':{str(k):nstr(v,125) for k,v in I.items()},
             'K2':{str(k):nstr(v,125) for k,v in K2.items()},
             'K3':{str(k):nstr(v,125) for k,v in K3.items()}}, open('ckpt_integrals135.pkl','wb'))
print("saved")
