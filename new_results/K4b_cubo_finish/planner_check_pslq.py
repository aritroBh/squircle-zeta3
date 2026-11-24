import json, mpmath as mp
mp.mp.dps=50
d=json.load(open('tetra_Sdk_final.json')); S=[mp.mpf(str(x)) for x in d['Sdk']]
r3=mp.sqrt(3); pi=mp.pi; z3=mp.zeta(3); L2=mp.log(2); L3=mp.log(3)
def series_mul(a,b,n=5): return [mp.fsum(a[i]*b[k-i] for i in range(k+1)) for k in range(n)]
def norm_series(c,n=5):
    # coefficients of exp(-c*log2* d/(1+d)) = exp(-c L2 (d - d^2 + d^3 - d^4))
    return mp.taylor(lambda x: mp.e**(-c*L2*x/(1+x)), 0, n-1)
for c in [1,2,3]:
    T=series_mul(S,norm_series(c))
    print(f'--- S(1+d) = 2^({c} d/(1+d)) T(d):  T0..T3 =',[mp.nstr(t,15) for t in T[:4]])
    for k,names,vals in [(1,['1','s3'],[1,r3]),(2,['1','s3','pi2','s3pi2'],[1,r3,pi**2,r3*pi**2]),
                         (3,['1','s3','pi2','s3pi2','z3','s3z3'],[1,r3,pi**2,r3*pi**2,z3,r3*z3])]:
        r=mp.pslq([T[k]]+vals,tol=mp.mpf(10)**-44,maxcoeff=10**5,maxsteps=10**6)
        if r: print(f'   T{k}:',' + '.join(f'{x}*{n}' for x,n in zip(r,['T']+names) if x),' res',mp.nstr(mp.fsum(x*v for x,v in zip(r,[T[k]]+vals)),3))
        else: print(f'   T{k}: none')
T=series_mul(S,norm_series(2))
names=['1','s3','pi2','s3pi2','z3','s3z3','pi4','s3pi4']; vals=[1,r3,pi**2,r3*pi**2,z3,r3*z3,pi**4,r3*pi**4]
r=mp.pslq([T[4]]+vals,tol=mp.mpf(10)**-46,maxcoeff=10**5,maxsteps=10**6)
print('T4:', (' + '.join(f'{x}*{n}' for x,n in zip(r,['T']+names) if x)+'  res '+mp.nstr(mp.fsum(x*v for x,v in zip(r,[T[4]]+vals)),3)) if r else 'none')
S1=(9+3*r3)/4
print('check T3 closed form:', mp.nstr(T[3]-(-23*S1/48+(7-r3)*pi**2/64-r3*z3/6),3), ' T2:', mp.nstr(T[2]-(-S1/8+r3*pi**2/96),3))
