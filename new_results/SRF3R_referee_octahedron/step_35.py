from mpmath import mp, mpf, log, exp, factorial
mp.dps = 40
prev = None; mono = True; worstV=None
for i in range(20, 20001, 37):
    V = mpf(i)/1000
    c = 2*exp(V)*abs(log(1-exp(-V)))
    if prev is not None and c > prev + mpf('1e-30'): mono=False; worstV=V
    prev = c
print("C(V) monotone decreasing on [0.02,20]:", mono)
for V in ['0.2','1','5','20']:
    V=mpf(V); print(f"  C({V}) = {mp.nstr(2*exp(V)*abs(log(1-exp(-V))),10)}")
ok1 = all(factorial(n-j)/factorial(n) <= 1/factorial(j) for n in range(1,60) for j in range(1,n+1))
ok2 = all((j+1)**(n-j+1) >= 2**(n-j+1) for n in range(1,60) for j in range(1,n+1))
print("(n-j)!/n! <= 1/j! (1<=j<=n<60):", ok1, "; (j+1)^{n-j+1} >= 2^{n-j+1}:", ok2)
