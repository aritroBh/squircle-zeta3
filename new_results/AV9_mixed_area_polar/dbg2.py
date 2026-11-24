
import sympy as sp, time
import av9core as C
lv, lw, t = C.lv, C.lw, C.t
s1, s2, v = sp.symbols('s1 s2 v')
N = 5
T=time.time()
def log(*a): print(f"[{time.time()-T:6.1f}]", *a, flush=True)
w = 2*(t*s1-t*s2)*sum((2*t*s1)**m for m in range(N+1))
w = sp.expand(w); log("w")
S1 = 0
for k in range(1, N+1):
    S1 += (v*lv**k + (1-v)*lw**k)/sp.factorial(k) * w**k
    S1 = sp.expand(S1)
    S1 = sum(S1.coeff(t,m)*t**m for m in range(1, N+1))
    log("S1 k=",k)
lnS = C.series_log1p_var(S1, t, N); log("lnS")
logp = sum(sp.Rational((-1)**(m+1),m)*(2*t*s2)**m for m in range(1,N+1))
bracket = logp - t*s2*(lv+lw) + (sp.Rational(1,2)-t*s1)*sum(lnS[k]*t**k for k in range(1,N+1))
bracket = sp.expand(bracket); log("bracket size", sp.count_ops(bracket))
bracket = sum(bracket.coeff(t,m)*t**m for m in range(1,N+1))
e = C.series_exp_var(bracket, t, N); log("exp")
for n in range(1,N+1):
    log(f"e[{n}] ops:", sp.count_ops(e[n]))
