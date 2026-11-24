
pslq2 = r'''
import mpmath as mp, sympy as sp, pickle, time, sys
c_exact = {int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_selfdual.pkl','rb')).items()}
def basis_elements(N):
    out = []
    for a in range(N+1):
        for b in range(N//2+1):
            for c in range(N//3+1):
                for dd in range(N//5+1):
                    for e in range(N//7+1):
                        if 1 <= a+2*b+3*c+5*dd+7*e <= N: out.append((a,b,c,dd,e))
    return out
def exact_vec(n, basis):
    Z2 = sp.Symbol('Z2'); dpi = sp.Symbol('dpi')
    e = sp.expand(sp.expand(c_exact[n]).subs(sp.pi, dpi).subs(dpi**2, 6*Z2))
    vec = []
    for m in basis:
        cf = sp.expand(e).coeff(Z2,m[1]).coeff(sp.log(2),m[0]).coeff(sp.zeta(3),m[2]).coeff(sp.zeta(5),m[3]).coeff(sp.zeta(7),m[4])
        vec.append(cf)
    return vec
n, dps = int(sys.argv[1]), int(sys.argv[2])
mp.mp.dps = dps
L2 = mp.log(2); z = {k: mp.zeta(k) for k in (2,3,5,7)}
target = mp.mpf(sp.N(c_exact[n], dps+10))
basis = basis_elements(n)
vals = [L2**m[0]*z[2]**m[1]*z[3]**m[2]*z[5]**m[3]*z[7]**m[4] for m in basis]
t0=time.time()
rel = mp.pslq([target]+vals, maxcoeff=10**8, maxsteps=20000)
out = open(f'scratchpad/pslq_c{n}_{dps}.txt','w')
if rel is None:
    out.write("FAILED"); print("FAILED"); sys.exit()
resid = abs(rel[0]*target + sum(r*v for r,v in zip(rel[1:],vals)))
vec = exact_vec(n, basis)
i0 = next(i for i,v in enumerate(vec) if v)
rat = sp.Rational(rel[1:][i0], vec[i0])
match = all(sp.Rational(r) == rat*v for r,v in zip(rel[1:],vec))
msg = f"c_{n} dps={dps}: resid={mp.nstr(resid,3)} max|c|={max(abs(r) for r in rel)} dim={len(basis)+1} match={match} time={time.time()-t0:.0f}s"
out.write(msg); out.close()
print(msg)
'''
open('pslq2.py','w').write(pslq2)
import subprocess
for n,d in [(6,200),(7,200)]:
    r = subprocess.run(['<project>/.venv/bin/python','pslq2.py',str(n),str(d)],capture_output=True,text=True,timeout=1500)
    print(r.stdout, r.stderr[-200:])
