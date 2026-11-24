import json, mpmath as mp
mp.mp.dps=50
S=[mp.mpf(x) for x in json.load(open('modelA_h96.json'))['Sdk']]
s5=mp.sqrt(5); L2=mp.log(2); L3=mp.log(3); pi=mp.pi; z3=mp.zeta(3)
S0=2*s5/3; print('S0 check', mp.nstr(S[0]-S0,3))
q=S[1]/S[0]; print('q_A - (13/15 + (13 ln3 - 6 ln2)/5) =', mp.nstr(q-(mp.mpf(13)/15+(13*L3-6*L2)/5),3))
lc=(13*L3-6*L2)/5   # log c_A
T=[x for x in mp.taylor(lambda x: mp.e**(-lc*x/(1+x)),0,4)]
U=[mp.fsum(S[i]*T[k-i] for i in range(k+1)) for k in range(5)]   # S * c^{-d/(1+d)}
print('T1/S0 =', mp.nstr(U[1]/S0,30))
def hunt(name,X,names,vals,tol,mc):
    r=mp.pslq([X]+vals,tol=mp.mpf(10)**-tol,maxcoeff=mc,maxsteps=10**6)
    if r: print(name,':',' + '.join(f'{c}*{n}' for c,n in zip(r,['X']+names) if c),' res',mp.nstr(mp.fsum(c*v for c,v in zip(r,[X]+vals)),3))
    else: print(name,': none up to',mc)
b2=['1','pi2','L2^2','L2L3','L3^2','L2','L3']; v2=[1,pi**2,L2**2,L2*L3,L3**2,L2,L3]
hunt('T2/S0',U[2]/S0,b2,v2,44,10**5)
hunt('T2/S0 (pi2 only)',U[2]/S0,['1','pi2'],[1,pi**2],44,10**7)
hunt('T2/(S0 s5)',U[2]/S0,['1','pi2','s5pi2','s5'],[1,pi**2,s5*pi**2,s5],44,10**6)
b3=['1','pi2','z3','pi2L2','pi2L3']; v3=[1,pi**2,z3,pi**2*L2,pi**2*L3]
hunt('T3/S0',U[3]/S0,b3,v3,42,10**5)
hunt('T3/S0 +s5',U[3]/S0,b3+['s5','s5pi2','s5z3'],v3+[s5,s5*pi**2,s5*z3],42,10**4)
