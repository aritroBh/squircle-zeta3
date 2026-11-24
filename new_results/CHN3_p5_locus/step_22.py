exec(open('common.py').read())
import time
C4, Q7 = get_C4_Q7()
C4_1 = sp.expand(C4.subs(a3,1)); Q7_1 = sp.expand(Q7.subs(a3,1))
R_test = sp.resultant(C4_1.subs(a4,2), Q7_1.subs(a4,2), a5)
t=time.time()
f = sp.factor_list(R_test, a1, a2)
print("scalar:", f[0])
for g,m in f[1]:
    P = sp.Poly(g, a1, a2)
    print(f"  total-degree {P.total_degree()}, mult {m}: {str(g)[:100]}")
print(f"{time.time()-t:.1f}s")
