import mpmath, pickle
from mpmath import mp, mpf, log, exp, quad, nstr
mp.dps = 90
g1 = lambda s: log(1+exp(-2*s))
g1p = lambda s: -2*exp(-2*s)/(1+exp(-2*s))
def Iab(a,b): return quad(lambda s: s**a*g1(s)**b, [0,mpf('0.5'),1,2,4,8,16,40])
def K2int(m): return quad(lambda s: s**m*g1(s)*g1p(s)**2, [0,mpf('0.5'),1,2,4,8,16,40])
ck = pickle.load(open('ckpt_integrals90.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ck['I'].items()}
K2 = {eval(k): mpf(v) for k,v in ck['K2'].items()}
for key in [(1,7),(3,5),(4,4),(5,3),(6,2)]:
    I[key] = Iab(*key); print(key, nstr(I[key],40))
K2[7] = K2int(7); print('K2(7)', nstr(K2[7],40))
pickle.dump({'I':{str(k):nstr(v,80) for k,v in I.items()},
             'K2':{str(k):nstr(v,80) for k,v in K2.items()},
             'K3':ck['K3']}, open('ckpt_integrals90.pkl','wb'))
print("saved")
