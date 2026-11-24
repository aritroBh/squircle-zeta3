code = r'''
# step_12.py : triangle-end coefficients: careful PSLQ with independent basis, two precisions.
import mpmath as mp

def compute_cs(dps, K=6):
    mp.mp.dps = dps
    def f(m): return mp.gamma(m+1)**2/mp.gamma(2*m+1)
    def gder(m, j): return 2*mp.psi(j-1, m+1) - mp.mpf(2)**j*mp.psi(j-1, 2*m+1)
    def bell(fderivs):
        k = len(fderivs)
        B = [[mp.mpf(0)]*(k+1) for _ in range(k+1)]
        B[0][0] = mp.mpf(1)
        for n in range(1,k+1):
            for j in range(1,n+1):
                B[n][j] = sum(mp.binomial(n-1,i-1)*fderivs[i-1]*B[n-i][j-1] for i in range(1,n-j+2))
        return [B[n][1] for n in range(1,k+1)]
    cs = [mp.mpf(0)]*(K+1)
    for m in range(1, 700):
        fm = f(m)
        if fm < mp.mpf(10)**(-dps-10) and m > 50: break
        gs = [gder(m,j) for j in range(1,K+1)]
        Bs = bell(gs)
        for k in range(0,K+1):
            term = mp.mpf(m)**(k-2)*(fm*Bs[k-1] if k>0 else fm)
            cs[k] += (-1)**k/mp.factorial(k)*term
    return cs

cs100 = compute_cs(100)
rt3 = mp.sqrt(3)
Cl = mp.clsin(2, mp.pi/3)   # Clausen Cl_2(pi/3) = 1.0149416064096536...
print('Cl2(pi/3) =', mp.nstr(Cl, 25))
pool = [mp.mpf(1), mp.pi, mp.log(2), mp.log(3), mp.zeta(2), mp.zeta(3), mp.pi/rt3, Cl]
pnames = ['1','pi','log2','log3','z2','z3','pi/rt3','Cl2(pi/3)']
for k in [1,2,3]:
    rel = mp.pslq([cs100[k]]+pool, tol=mp.mpf('1e-40'), maxcoeff=10**8, maxsteps=10000)
    print(f'c_{k} PSLQ vs pool:', rel)
    rel = mp.pslq([cs100[k]]+pool+[Cl*mp.pi/rt3, mp.zeta(2)*mp.log(2)], tol=mp.mpf('1e-40'), maxcoeff=10**8, maxsteps=10000)
    print(f'c_{k} PSLQ vs pool+prods:', rel)

# also try z3-free statement for c1,c2: PSLQ against {1,pi,log2,log3,z2,pi/rt3,Cl} (no z3)
pool7 = [mp.mpf(1), mp.pi, mp.log(2), mp.log(3), mp.zeta(2), mp.pi/rt3, Cl]
for k in [1,2]:
    rel = mp.pslq([cs100[k]]+pool7, tol=mp.mpf('1e-40'), maxcoeff=10**8, maxsteps=10000)
    print(f'c_{k} PSLQ (no-z3 pool):', rel)
'''
open('steps/step_12.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_12.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-800:] if r.returncode else '')
open('progress.log','a').write('step_12: triangle-end PSLQ with Clausen pool\n')
