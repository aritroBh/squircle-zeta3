import mpmath as mp
mp.mp.dps = 100
print("dps:", mp.mp.dps)
z2=mp.zeta(2); z3=mp.zeta(3); z5=mp.zeta(5); pi=mp.pi

# Exact Beta moments via series: int_0^1 (log t)^a (log(1-t))^b dt
# (a,b)=(2,1): -2 * sum 1/(n(n+1)^3) -> we derived = -6+2z2+2z3. Verify:
I21 = -6 + 2*z2 + 2*z3
I21q = mp.quad(lambda t: mp.log(t)**2*mp.log(1-t), [0,1], method='tanh-sinh')
print("(2,1) exact vs quad resid:", mp.nstr(abs(I21-I21q),3))

# (a,b)=(4,1): -24 * sum_{n>=1} 1/(n(n+1)^5)
S41 = mp.nsum(lambda n: 1/(n*(n+1)**5), [1, mp.inf])
I41_series = -24*S41
I41q = mp.quad(lambda t: mp.log(t)**4*mp.log(1-t), [0,1], method='tanh-sinh')
print("(4,1) series vs quad resid:", mp.nstr(abs(I41_series-I41q),3))

# (a,b)=(3,2): -12 * sum_{n>=2} H_{n-1}/(n(n+1)^4)
def H(n): return mp.fsum([mp.mpf(1)/k for k in range(1,n+1)])
S32 = mp.nsum(lambda n: H(int(n)-1)/(n*(n+1)**4), [2, mp.inf])
I32_series = -12*S32
I32q = mp.quad(lambda t: mp.log(t)**3*mp.log(1-t)**2, [0,1], method='tanh-sinh')
print("(3,2) series vs quad resid:", mp.nstr(abs(I32_series-I32q),3))

# PSLQ against weight<=5 basis {1, z2, z3, z2^2, z2*z3, z5}
basis = [mp.mpf(1), z2, z3, z2**2, z2*z3, z5]
bn = ['1','z2','z3','z2^2','z2*z3','z5']
for name, I in [('(2,1)',I21q),('(4,1)',I41q),('(3,2)',I32q)]:
    rel = mp.pslq(mp.matrix([I]+basis), tol=mp.mpf('1e-90'), maxcoeff=10**8)
    print(f"\nBeta {name}: value={mp.nstr(I,20)}")
    if rel is not None and rel[0]!=0:
        rel=list(rel)
        expr = {bn[k]: -mp.mpf(rel[k+1])/rel[0] for k in range(len(basis))}
        print("  PSLQ: I =", {k: mp.nstr(v,8) for k,v in expr.items() if v!=0})
        print("  z5 coeff:", mp.nstr(-mp.mpf(rel[6])/rel[0],10), " z2z3 coeff:", mp.nstr(-mp.mpf(rel[5])/rel[0],10))
    else:
        print("  PSLQ failed:", rel)
