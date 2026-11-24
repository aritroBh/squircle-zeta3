import pickle, mpmath as mp, time
from q1s_exact import femb
t0=time.time()
mp.mp.dps=1150
r2=mp.sqrt(2); pi=mp.pi
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
R=pickle.load(open('Rreal2_V3p_1250.pkl','rb'))
Rre=mp.mpf(R['Rre'])
tests=[ [('z3',z3)], [('Lchi8',Lchi8)], [('pi3',pi**3)],
        [('z3',z3),('Lchi8',Lchi8)], [('z3',z3),('pi3',pi**3)], [('Lchi8',Lchi8),('pi3',pi**3)],
        [('z3',z3),('Lchi8',Lchi8),('pi3',pi**3)] ]
for basis in tests:
    vec=[Rre]+[v for _,v in basis]
    rel=mp.pslq(vec, tol=mp.mpf('1e-1090'), maxcoeff=10**60, maxsteps=3000)
    ok = rel and rel[0]!=0
    print([nm for nm,_ in basis], "->", rel if ok else "(none)")
print("time %.0fs"%(time.time()-t0))
