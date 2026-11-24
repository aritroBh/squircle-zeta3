import mpmath, pickle, time
from mpmath import mp, mpf, log, exp, quad, nstr
mp.dps = 210
g1 = lambda s: log(1+exp(-2*s))
g1p = lambda s: -2*exp(-2*s)/(1+exp(-2*s))
PTS = [0, mpf('0.5'), 1, 2, 4, 8, 16, 32, 64, 100, 125]
def Iab(a,b): return quad(lambda s: s**a*g1(s)**b, PTS)
def K2int(m): return quad(lambda s: s**m*g1(s)*g1p(s)**2, PTS)
need = [(4,3),(2,5),(0,7),(5,2),(1,6),(3,4),(0,8),(2,6),(6,2),(5,3),(4,4),(3,5),(1,7)]
t0=time.time()
I = {key: Iab(*key) for key in need}
K2 = {6: K2int(6), 7: K2int(7)}
pickle.dump({'I':{str(k):nstr(v,200) for k,v in I.items()},'K2':{str(k):nstr(v,200) for k,v in K2.items()}},
            open('ckpt_integrals210.pkl','wb'))
print("time", time.time()-t0, "s")
print("sample I(5,2) =", nstr(I[(5,2)],50))
