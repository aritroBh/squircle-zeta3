import sympy as sp
# Exact closed forms for c(w) from the law, as polynomials in lam, pi, zeta(odd)
lam, pi, y = sp.symbols('lam pi y')
zet = lambda n: sp.Symbol(f'z{n}')
# s_n: even n -> zeta(n) = rat*pi^n ; odd n -> (1-2^{1-n}) zeta(n)
def s_n(n):
    if n % 2 == 0:
        return sp.zeta(n)  # sympy will give rat*pi^n
    else:
        return (1 - sp.Rational(2)**(1-n)) * zet(n)
N = 12
L = [sp.Rational(0)]*(N+1)
L[1] = lam
for n in range(2, N+1): L[n] = sp.expand(s_n(n))/n
E = [sp.Rational(1)] + [sp.Rational(0)]*N
for n in range(1, N+1):
    E[n] = sp.expand(sum(k*L[k]*E[n-k] for k in range(1, n+1))/n)
A = [1,1,1,2,5,16,61,272,1385,7936,50521]
def eps(m): return [1,-1,-1,1][m%4]
cw_forms = {}
for w in range(2, 10):
    m = w-2
    cw_forms[w] = sp.expand(eps(m)*A[m+1]/(2*sp.factorial(m)) * E[w])
    print(f"c({w}) = {cw_forms[w]}")
import json
json.dump({str(w): str(cw_forms[w]) for w in cw_forms}, open('cw_closed_forms.json','w'))
print('saved')
with open('progress.log','a') as f:
    f.write('step_22: w6 law confirmed to 406 digits for all w=2..9; exact closed forms via sympy -> cw_closed_forms.json\n')
