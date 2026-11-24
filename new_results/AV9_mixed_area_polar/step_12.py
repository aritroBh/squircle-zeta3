
# step_12.py : PSLQ record for c_5..c_8 at dps 100 and 150; incomplete-basis control
import mpmath as mp
import sympy as sp
import pickle

c_exact = {int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_selfdual.pkl','rb')).items()}

def basis_elements(N):
    # monomials ln2^a zeta(2)^b zeta(3)^c zeta(5)^d zeta(7)^e with weight <= N, plus 1
    out = []
    for a in range(N+1):
        for b in range(N//2+1):
            for c in range(N//3+1):
                for d in range(N//5+1):
                    for e in range(N//7+1):
                        w = a+2*b+3*c+5*d+7*e
                        if w <= N:
                            out.append((a,b,c,d,e))
    return out

def run_pslq(n, dps, drop=None):
    mp.mp.dps = dps
    L2 = mp.log(2); z = {k: mp.zeta(k) for k in (2,3,5,7)}
    target = mp.mpf(sp.N(c_exact[n], dps+10))
    basis = basis_elements(n)
    if drop: basis = [m for m in basis if not any(m[i]==drop[i] for i in range(5) if drop[i])]
    vals = []
    for (a,b,c,d,e) in basis:
        vals.append(L2**a * z[2]**b * z[3]**c * z[5]**d * z[7]**e)
    rel = mp.pslq([target]+vals, maxcoeff=10**7, maxsteps=10000)
    return basis, rel

for n in [5,6,7,8]:
    for dps in [100, 150]:
        basis, rel = run_pslq(n, dps)
        if rel is None:
            print(f"c_{n} dps={dps}: PSLQ FAILED")
            continue
        # rebuild candidate expression
        expr = -rel[0]**-1 * sum(sp.Integer(r)*sp.log(2)**m[0]*sp.zeta(2)**m[1]*sp.zeta(3)**m[2]*sp.zeta(5)**m[3]*sp.zeta(7)**m[4] for r,m in zip(rel[1:],basis))
        match = sp.simplify(expr - c_exact[n]) == 0
        print(f"c_{n} dps={dps}: max|coeff|={max(abs(r) for r in rel)}  exact-match={match}")
print("\n-- incomplete basis controls (dps=100) --")
for n, drop in [(5, (0,0,1,0,0)), (6, (0,0,0,1,0)), (8, (0,0,0,0,1))]:
    basis, rel = run_pslq(n, 100, drop=drop)
    tag = { (0,0,1,0,0):'no zeta(3)', (0,0,0,1,0):'no zeta(5)', (0,0,0,0,1):'no zeta(7)'}[drop]
    print(f"c_{n} {tag}: PSLQ result = {'FAILED (no relation <1e7)' if rel is None else 'FOUND maxcoeff '+str(max(abs(r) for r in rel))}")
