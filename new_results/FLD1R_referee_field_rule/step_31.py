import sys, json, time, pickle
sys.path.insert(0,'.')
from fractions import Fraction as Fr
sys.argv=['x','7']
src=open('step_22_musymbolic.py').read(); ns={}
exec(src[:src.index("if __name__")], ns)
KMAX=ns['KMAX']; mu_formula=ns['mu_formula']
print('KMAX =', KMAX)
t0=time.time()
F7 = mu_formula(7)
print('mu_formula(7) done in', round(time.time()-t0,1),'s;', len(F7.terms()), 'monomials', flush=True)
pickle.dump([(mon, str(cf.as_expr())) for mon,cf in F7.terms()], open('mu_symbolic_K7.pkl','wb'))
# verify against exact table n<=30
MU = json.load(open('mu_exact_K7_N30.json'))
def Hf(r,m): return sum((Fr(1,i**r) for i in range(1,m+1)),Fr(0))
def evalK(c,nv): return Fr(int(c.numer(nv)),1)/Fr(int(c.denom(nv)),1)
ok=True
for nv in range(1,31):
    acc={}
    for mon,cf in F7.terms():
        cval=evalK(cf,nv); hv=Fr(1)
        for r in range(1,KMAX+1):
            hv*=Hf(r,nv-1)**mon[r-1]; hv*=Hf(r,2*nv-1)**mon[KMAX+r-1]
        ze=tuple(list(mon[2*KMAX:])+[0]*(6-len(mon[2*KMAX:])))[:6]
        zk=','.join(map(str,ze)); acc[zk]=acc.get(zk,Fr(0))+cval*hv
    ref=MU['7'][str(nv)]
    for kk in set(list(acc)+list(ref)):
        if Fr(acc.get(kk,'0'))!=Fr(ref.get(kk,'0')): ok=False; print('MISMATCH',nv,kk)
print('mu_formula(7) == exact table (n<=30):', ok)
open('progress.log','a').write(f'step04a: symbolic-n mu_formula(7) ({len(F7.terms())} monomials) == exact table n<=30: {ok}\n')
