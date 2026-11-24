import sys, time, pickle
sys.path.insert(0,'.')
import zsum
from zsum import K, nn, ONE, add_into
sys.argv=['x','7']
src=open('step_22_musymbolic.py').read(); ns={}
exec(src[:src.index("if __name__")], ns)
KMAX=ns['KMAX']; mu_formula=ns['mu_formula']
t0=time.time()
F7 = mu_formula(7)
def harm_n_minus_1(r): return {(r,): K(2**(r-1)), (-r,): K(2**(r-1)), (): -ONE/nn**r}
def harm_2n_minus_1(r): return {(r,): ONE, (): -ONE/(2*nn)**r}
def to_zsum(F):
    parts={}
    for mon,cf in F.terms():
        E={(): ONE}
        for r in range(1,KMAX+1):
            for _ in range(mon[r-1]): E=zsum.mul(E,harm_n_minus_1(r))
            for _ in range(mon[KMAX+r-1]): E=zsum.mul(E,harm_2n_minus_1(r))
        zkey=tuple(list(mon[2*KMAX:])+[0]*(6-len(mon[2*KMAX:])))[:6]
        parts.setdefault(zkey,{})
        for w,c in E.items(): add_into(parts[zkey],w,cf*c)
    return parts
Z=to_zsum(F7)
print('to_zsum(7) in', round(time.time()-t0,1),'s; words per zeta-monomial:', {zk:len(E) for zk,E in sorted(Z.items())}, flush=True)
pickle.dump({zk:[(w,c.as_expr()) for w,c in E.items()] for zk,E in Z.items()}, open('mu_zsum_K7.pkl','wb'))
print('saved mu_zsum_K7.pkl')
open('progress.log','a').write(f'step04b: Z-sum canonical form k=7: {sum(len(E) for E in Z.values())} word-terms over {len(Z)} zeta-monomials; saved mu_zsum_K7.pkl\n')
