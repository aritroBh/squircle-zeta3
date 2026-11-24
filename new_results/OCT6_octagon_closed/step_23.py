import pickle, sympy as sp, mpmath as mp
import mixed_engine as ME
mp.mp.dps = 750
r2m=mp.sqrt(2); k1m=r2m-1; k2m=r2m/2-1
x1m=-k1m; x2m=-k2m
def AA_of(xv):
    # sum H_n xv^n / n^2 ; |xv|<1 ; need ~ 750/log10(1/|xv|) terms
    import math
    N = int(760/mp.log10(1/abs(xv)))+50
    s=mp.mpf(0); Hn=mp.mpf(0); xp=mp.mpf(1)
    for nn in range(1,N+1):
        Hn+=1/mp.mpf(nn); xp*=xv; s+=xp*Hn/nn**2
    return s
def Li21(a,b):
    # sum_{n>m>=1} a^n b^m/(n^2 m); |a|<1
    N = int(760/mp.log10(1/abs(a)))+50
    s=mp.mpf(0); ap=mp.mpf(1); inner=mp.mpf(0); bp=mp.mpf(1)
    for n in range(2,N+1):
        ap*=a; bp*=b; inner += bp/(n-1)
        s += ap*inner/n**2
    return s
import time; t0=time.time()
num = {
 ME.l2: mp.log(2), ME.lu: mp.log(1+r2m), ME.z2: mp.pi**2/6, ME.z3: mp.zeta(3),
 ME.L2_1: mp.polylog(2,x1m), ME.L3_1: mp.polylog(3,x1m),
 ME.L2_2: mp.polylog(2,x2m), ME.L3_2: mp.polylog(3,x2m),
 ME.L2k: mp.polylog(2,k1m), ME.L3k: mp.polylog(3,k1m),
 ME.L2d: mp.polylog(2,2-r2m), ME.L3d: mp.polylog(3,2-r2m),
}
print("polylogs done", time.time()-t0); t0=time.time()
num[ME.AA_1]=AA_of(x1m); num[ME.AA_2]=AA_of(x2m)
num[ME.Lam]=Li21(x1m,-r2m/2); num[ME.swap]=Li21(x2m,-r2m)
print("series done", time.time()-t0)
# cross-check vs atoms600 (600-digit values)
atoms = {k: mp.mpf(v) for k,v in pickle.load(open('oct5/atoms600.pkl','rb')).items()}
for sym, key in [(ME.AA_1,'AA_1'),(ME.AA_2,'AA_2'),(ME.Lam,'Li21_a'),(ME.swap,'Li21_b'),(ME.L2_1,'L2_1'),(ME.L3_1,'L3_1')]:
    print(key, "diff vs atoms600:", mp.nstr(abs(num[sym]-atoms[key]),3))
pickle.dump({str(k): mp.nstr(v, 745) for k,v in num.items()}, open('num745.pkl','wb'))
print("saved num745.pkl")
