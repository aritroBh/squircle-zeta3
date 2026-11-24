
import os, sys, time
os.chdir('<project>/research/08_theory/CHN3_p5_locus')
sys.path.insert(0, '.')
exec(open('common.py').read())
import sympy as sp
t0=time.time()
C4, Q7 = get_C4_Q7()
C4_1 = sp.expand(C4.subs(a3,1)); Q7_1 = sp.expand(Q7.subs(a3,1))
print(f"setup {time.time()-t0:.1f}s; starting resultant", flush=True)
t0=time.time()
R = sp.resultant(C4_1, Q7_1, a5)
print(f"resultant done {time.time()-t0:.1f}s", flush=True)
R = sp.expand(R)
open('R_a1a2a4.txt','w').write(str(R))
P = sp.Poly(R, a1,a2,a4)
print("R total degree:", P.total_degree(), " terms:", len(P.terms()), flush=True)
t0=time.time()
f = sp.factor_list(R, a1,a2,a4)
print(f"factor done {time.time()-t0:.1f}s; scalar {f[0]}", flush=True)
with open('R_factors.txt','w') as fh:
    fh.write(f"scalar: {f[0]}\n")
    for g,m in f[1]:
        Pg=sp.Poly(g,a1,a2,a4)
        fh.write(f"mult {m} deg {Pg.total_degree()} terms {len(Pg.terms())}: {g}\n")
        print(f"  mult {m} deg {Pg.total_degree()} terms {len(Pg.terms())}", flush=True)
print("DONE", flush=True)
