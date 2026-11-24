import mpmath as mp
mp.mp.dps = 60

def Li21(a, c, N=500):
    a, c = mp.mpf(a), mp.mpf(c)
    return mp.fsum([a**n * c**m / (mp.mpf(n)**2 * m) for n in range(2, N) for m in range(1, n)])

def reduce_or_survive(target, basis_vals, basis_names, dps=None, tol_exp=None, maxcoeff=10**10):
    """Iteratively PSLQ. Internal relations (not involving target) are quotiented
    (one element dropped). Returns ('REDUCES', rel, names) or ('SURVIVES', None, dropped_internal)."""
    vals = list(basis_vals); names = list(basis_names)
    internal = []
    while True:
        v = vals + [target]
        nm = names + ['TARGET']
        tol = mp.mpf(10)**(-(mp.mp.dps-8)) if tol_exp is None else mp.mpf(10)**(-tol_exp)
        rel = mp.pslq(v, tol=tol, maxcoeff=maxcoeff, maxsteps=2000)
        if rel is None:
            return ('SURVIVES', None, internal, len(v))
        tpos = len(v)-1
        if rel[tpos] != 0:
            return ('REDUCES', rel, nm, internal)
        # internal relation: drop the last basis element with nonzero coefficient
        nz = [i for i in range(tpos) if rel[i] != 0]
        drop = nz[-1]
        internal.append((names[drop], rel, list(names)))
        vals.pop(drop); names.pop(drop)

x1 = mp.mpf(-1)/3; x2 = mp.mpf(-1)/9
b1 = mp.mpf(1)/4; b2 = mp.mpf(1)/10
Lam = [Li21(x1, x2/x1), Li21(x2, x1/x2), Li21(b1, b2/b1), Li21(b2, b1/b2)]
l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
z2 = mp.pi**2/6; z3 = mp.zeta(3)
def Li2(v): return mp.polylog(2, mp.mpf(v))
def Li3(v): return mp.polylog(3, mp.mpf(v))
# independent-ish basis (avoid Landen pairs): logs 2,3,5; Li2/Li3 at x1,x2,b1,b2,-1/2,-2/3,2/3,3/5
pts = [mp.mpf(-1)/3, mp.mpf(-1)/9, mp.mpf(1)/4, mp.mpf(1)/10, mp.mpf(-1)/2, mp.mpf(-2)/3, mp.mpf(2)/3, mp.mpf(3)/5]
pnames = ["-1/3","-1/9","1/4","1/10","-1/2","-2/3","2/3","3/5"]
basis_vals = [mp.mpf(1), l2, l3, l5, z2]
basis_names = ["1","log2","log3","log5","z2"]
for p,nm in zip(pts,pnames):
    basis_vals.append(Li2(p)); basis_names.append(f"Li2({nm})")
for a,na in zip([l2,l3,l5],["log2","log3","log5"]):
    for b,nb in zip([l2,l3,l5],["log2","log3","log5"]):
        basis_vals.append(a*b); basis_names.append(f"{na}{nb}")
basis_vals.append(z3); basis_names.append("z3")
for p,nm in zip(pts,pnames):
    basis_vals.append(Li3(p)); basis_names.append(f"Li3({nm})")
for a,na in zip([l2,l3,l5],["log2","log3","log5"]):
    basis_vals.append(z2*a); basis_names.append(f"z2{na}")
    for p,nm2 in zip(pts,pnames):
        basis_vals.append(Li2(p)*a); basis_names.append(f"Li2({nm2}){na}")
import itertools
for a,na in zip([l2,l3,l5],["log2","log3","log5"]):
    basis_vals.append(a**3); basis_names.append(f"{na}^3")
for (a,na),(b,nb) in itertools.combinations(list(zip([l2,l3,l5],["log2","log3","log5"])),2):
    basis_vals.append(a*a*b); basis_vals.append(a*b*b)
    basis_names.append(f"{na}2{nb}"); basis_names.append(f"{na}{nb}2")
print("basis size:", len(basis_vals))
names4 = ["Lam1=Li21(-1/3,1/3)","Lam2=Li21(-1/9,3)","Lam3=Li21(1/4,2/5)","Lam4=Li21(1/10,5/2)"]
for nm4, L in zip(names4, Lam):
    st, rel, nm, extra = reduce_or_survive(L, basis_vals, basis_names, maxcoeff=10**9)
    if st=='REDUCES':
        terms = [(nm[i], rel[i]) for i in range(len(rel)) if rel[i]!=0]
        print(f"{nm4}: REDUCES at 60dps:")
        for n_, c_ in terms: print(f"     {c_} * {n_}")
    else:
        print(f"{nm4}: SURVIVES at 60dps (basis {extra} after {len(rel) if rel else 0} internal quotients)" if False else f"{nm4}: SURVIVES at 60dps, {len(extra[2])} internal relations quotiented, final basis {extra[3]}")
