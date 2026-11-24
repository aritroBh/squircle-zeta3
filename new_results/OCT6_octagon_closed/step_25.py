import pickle, mpmath as mp, sympy as sp
import mixed_engine as ME
mp.mp.dps = 750
r2m=mp.sqrt(2); k1m=r2m-1; k2m=r2m/2-1
x1m=-k1m; x2m=-k2m
def Li21(a,b):
    N = int(760/mp.log10(1/abs(a)))+50
    s=mp.mpf(0); ap=a*a; inner=b; bp=b
    for n in range(2,N+1):
        if n>2: ap*=a; bp*=b; inner+=bp/(n-1)
        s += ap*inner/n**2
    return s
num = {sp.symbols(k): mp.mpf(v) for k,v in pickle.load(open('num745.pkl','rb')).items()}
Lam = Li21(x1m,-r2m/2); swap = Li21(x2m,-r2m)
atoms = {k: mp.mpf(v) for k,v in pickle.load(open('oct5/atoms600.pkl','rb')).items()}
print("Lam diff vs atoms600:", mp.nstr(abs(Lam-atoms['Li21_a']),3))
print("swap diff vs atoms600:", mp.nstr(abs(swap-atoms['Li21_b']),3))
num[ME.Lam]=Lam; num[ME.swap]=swap
pickle.dump({str(k): mp.nstr(v, 745) for k,v in num.items()}, open('num745.pkl','wb'))
print("num745.pkl fixed")
