exec(open('common.py').read())
import time
C4, Q7 = get_C4_Q7()
C4_1 = sp.expand(C4.subs(a3,1))
Q7_1 = sp.expand(Q7.subs(a3,1))
# Specialize a4=2 to gauge cost (univariate resultant in a5 over Q(a1,a2))
t=time.time()
R_test = sp.resultant(C4_1.subs(a4,2), Q7_1.subs(a4,2), a5)
print("deg in (a1,a2):", sp.Poly(R_test, a1,a2).total_degree(), "terms:", len(sp.Poly(R_test,a1,a2).terms()), f"  {time.time()-t:.1f}s")
t=time.time()
f = sp.factor_list(R_test)
print("factors:", [(sp.degree(g), m) for g,m in f[1]], f"  {time.time()-t:.1f}s")
for g,m in f[1]:
    print("  deg", sp.Poly(g,a1,a2).total_degree(), ":", str(g)[:120], " x", m)
