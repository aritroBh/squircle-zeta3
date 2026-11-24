import sympy as sp, time
t0 = time.time()
x = sp.Symbol('x')

# Stirling S(n,k) table
maxn = 25
S = [[0]*(maxn+1) for _ in range(maxn+1)]
S[0][0] = 1
for n in range(1, maxn+1):
    for k in range(1, n+1):
        S[n][k] = k*S[n-1][k] + S[n-1][k-1]

# P_5(x) = S(x, x-5): interpolate degree-10 polynomial through (n, S(n,n-5)), n=0..10
pts = [(n, S[n][n-5]) for n in range(0, 11)]
P5 = sp.interpolate(pts, x).expand()
print("P5 interpolated:", sp.factor(P5))
# leading coefficient should be 1/(2^5*5!) = 1/3840
print("leading:", sp.LC(sp.Poly(P5, x)), "1/3840 =", sp.Rational(1,3840))
# Verify against more values
ok_vals = all(P5.subs(x, n) == S[n][n-5] for n in range(0, 21))
print("matches S(n,n-5) for n=0..20:", ok_vals)

# P_4 from prove_S34.py
P4 = x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760
# PROOF 1: recurrence P5(x+1) = (x+1-5) P4(x) + P5(x)  -- exact symbolic identity
rec = sp.expand(P5.subs(x, x+1) - ((x+1-5)*P4 + P5))
print("recurrence P5(x+1) - ((x-4)P4+P5) == 0 :", rec == 0)
# PROOF 2: anchors P5(a) = 0 for a = 0..5
anchors = all(P5.subs(x, a) == 0 for a in range(0, 6))
print("anchors P5(a)=0, a=0..5:", anchors)
# Factored form for the proof document
print("P5 factored:", sp.factor(P5*sp.Rational(40320,1)))
print("time", round(time.time()-t0,1))
