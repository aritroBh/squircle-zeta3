import mpmath, pickle
from mpmath import mp, mpf, zeta, log, nstr
mp.dps = 120
F = {k: mpf(v) for k,v in pickle.load(open('ckpt_Fk_125.pkl','rb')).items()}
Z = {n: zeta(n) for n in range(2,12)}
# zeta-only monomial bases
bases = {
 4: [('Z2^2', Z[2]**2)],
 5: [('Z2Z3', Z[2]*Z[3]), ('Z5', Z[5])],
 6: [('Z2^3', Z[2]**3), ('Z3^2', Z[3]**2)],
 7: [('Z2^2Z3', Z[2]**2*Z[3]), ('Z2Z5', Z[2]*Z[5]), ('Z7', Z[7])],
 8: [('Z2^4', Z[2]**4), ('Z2Z3^2', Z[2]*Z[3]**2), ('Z3Z5', Z[3]*Z[5]), ('Z8', zeta(8))],
 9: [('Z2^3Z3', Z[2]**3*Z[3]), ('Z2^2Z5', Z[2]**2*Z[5]), ('Z2Z7', Z[2]*Z[7]),
     ('Z3^3', Z[3]**3), ('Z3Z5?', Z[3]*Z[5]*0), ('Z9', Z[9])],
}
bases[9] = bases[9][:4] + bases[9][5:]
for k in range(3,9):
    w = k+1
    B = bases[w]
    rel = mpmath.pslq([v for _,v in B]+[F[k]], maxcoeff=10**8)
    if rel and rel[-1]!=0:
        mc = max(abs(x) for x in rel); N = len(B)+1
        bound = mpf(10)**(mpf(115)/(N-1))
        txt = " + ".join(f"({r})*{nm}" for r,(nm,_) in zip(rel,B) if r) + f" + ({rel[-1]})*F{k}"
        print(f"w{w} N{N} maxc {nstr(mc,4)} bound {nstr(bound,3)}: {txt}")
    else:
        print(f"w{w}: NONE")
