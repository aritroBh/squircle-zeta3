import json, mpmath as mp
mp.mp.dps=70
S=[mp.mpf(x) for x in json.load(open('modelA_h192.json'))['Sdk']]
S128=[mp.mpf(x) for x in json.load(open('modelA_h128.json'))['Sdk']]
s5=mp.sqrt(5); L2=mp.log(2); L3=mp.log(3); pi=mp.pi; z3=mp.zeta(3); L=mp.log(mp.mpf(4)/3)
print('grid agreement h128 vs h192:', [mp.nstr(abs(a-b),3) for a,b in zip(S,S128)])
lc=(13*L3-6*L2)/5
T=mp.taylor(lambda x: mp.e**(-lc*x/(1+x)),0,4)
U=[mp.fsum(S[i]*T[k-i] for i in range(k+1)) for k in range(5)]
d2=mp.polylog(2,mp.mpf(-1)/3)
T2c=s5*(mp.mpf(-254)/675+mp.mpf(79)/8100*pi**2-mp.mpf(118)/75*d2-mp.mpf(4)/15*L+mp.mpf(76)/75*L**2)
print('T2 - claimed =', mp.nstr(U[2]-T2c,3))
X=U[3]/s5
names=['1','pi2','z3','Li3(-1/3)','Li3(1/4)','Li3(3/4)','Li2(-1/3)L','pi2L','L','L2','L3','Li2(-1/3)']
vals=[1,pi**2,z3,mp.polylog(3,mp.mpf(-1)/3),mp.polylog(3,mp.mpf(1)/4),mp.polylog(3,mp.mpf(3)/4),d2*L,pi**2*L,L,L**2,L**3,d2]
for tol in (52,50):
    r=mp.pslq([X]+vals,tol=mp.mpf(10)**-tol,maxcoeff=10**6,maxsteps=10**6)
    if r: print(f'T3/s5 (tol 1e-{tol}):',' + '.join(f'{c}*{n}' for c,n in zip(r,['X']+names) if c),' res',mp.nstr(mp.fsum(c*v for c,v in zip(r,[X]+vals)),3))
    else: print('T3: none')
