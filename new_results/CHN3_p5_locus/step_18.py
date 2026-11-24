exec(open('common.py').read())
C4, Q7 = get_C4_Q7()
import time
t=time.time()
fC4 = sp.factor_list(C4, a1,a2,a3,a4,a5)
print("factor C4:", [(str(g)[:60], m) for g,m in fC4[1]], " coeff:", fC4[0], f"{time.time()-t:.1f}s")
t=time.time()
fQ7 = sp.factor_list(Q7, a1,a2,a3,a4,a5)
print("factor Q7':", [(str(g)[:60], m) for g,m in fQ7[1]], " coeff:", fQ7[0], f"{time.time()-t:.1f}s")
# reversal symmetry
rev = {a1:a5, a2:a4, a3:a3, a4:a2, a5:a1}
print("C4 reversal-invariant:", sp.expand(C4 - C4.subs(rev))==0)
print("Q7' reversal-invariant:", sp.expand(Q7 - Q7.subs(rev))==0)
