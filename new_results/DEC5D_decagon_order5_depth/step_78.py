
import pickle
import mpmath as mp
mp.mp.dps = 250
r5 = mp.sqrt(5); phi = (1+r5)/2
z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
# exact value table: basis vectors per weight
ZV = [mp.mpf(1), z2, z3, z2**2, z5, z2*z3]
ZN = ['1','z2','z3','z2^2','z5','z2z3']
def recon(val, wmax):
    # value = sum c_k ZV[k] with weight(ZN[k]) <= wmax
    cand = [(ZN[k], ZV[k]) for k in range(6) if [0,2,3,4,5,5][k] <= wmax]
    vec = [val] + [c[1] for c in cand]
    rel = mp.pslq(vec, tol=mp.mpf('1e-200'), maxcoeff=10**25, maxsteps=10000)
    if rel is None: return None
    out = {}
    for i, (n, _) in enumerate(cand):
        if rel[i+1]: out[n] = -mp.mpf(rel[i+1])/mp.mpf(rel[0])
    # verify
    s = sum(out.get(ZN[k], 0)*ZV[k] for k in range(6) if ZN[k] in out)
    if abs(s - val) > mp.mpf('1e-200'): return None
    return out

parsed = pickle.load(open("parsed_mons_dec.pkl","rb"))
localkeys = set()
for mon, cu in parsed:
    a,b,c,d,e,f = mon
    if c==0 and d==0:
        for n in cu: localkeys.add((n,a,b))
exactI = {}
fails = []
for k in sorted(localkeys):
    n,a,b = k
    f = lambda x: x**n * mp.log(x)**a * mp.log(1-x)**b
    v = mp.quad(f, [0,1], method='tanh-sinh')
    r = recon(v, a+b)
    if r is None: fails.append(k)
    else: exactI[k] = r
print("reconstructed:", len(exactI), "failed:", fails)
for k in sorted(exactI)[:6]:
    print("  ", k, exactI[k])
pickle.dump({k: {n: str(c) for n, c in v.items()} for k, v in exactI.items()}, open("localI_exact.pkl","wb"))
