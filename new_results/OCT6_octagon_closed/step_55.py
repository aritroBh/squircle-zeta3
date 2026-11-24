import sympy as sp
x = sp.symbols('x')
# R1 as function identity on (-1,0): substitutions
# l2 = 2*log(1-x), lu = -log(-x), x2 = x/(x-1), r2h = 1/(1-x), d = 1+x, mr2h = 1/(x-1)
L2f = lambda z: sp.polylog(2, z)
L3f = lambda z: sp.polylog(3, z)
l2f = 2*sp.log(1-x); luf = -sp.log(-x)
x2f = x/(x-1); r2hf = 1/(1-x); df = 1+x
R1f = (sp.Rational(1,2)*sp.zeta(3) + sp.Rational(3,2)*L3f(x) + sp.Rational(1,2)*L3f(x2f) - sp.Rational(1,2)*L3f(r2hf)
       - sp.Rational(1,2)*sp.pi**2/6*l2f - sp.Rational(1,2)*l2f**2*luf + sp.Rational(1,2)*l2f*luf**2
       + sp.Rational(1,2)*L2f(x)*l2f - sp.Rational(1,2)*L2f(x)*luf - sp.Rational(1,2)*L2f(x2f)*luf + L2f(df)*l2f)
# A(x) = sum H_n x^n/n^2 ; A'(x) = (Li2(x) + log(1-x)^2/2)/x
Aprime = (L2f(x) + sp.log(1-x)**2/2)/x
F1 = sp.diff(R1f, x) - Aprime
F1e = sp.expand(F1)
print("F1'(x) expanded terms:", len(sp.Add.make_args(F1e)))
# Now collapse using FUNCTION identities for dilogs:
# (i) Landen: Li2(x/(x-1)) = -Li2(x) - log(1-x)^2/2
# (ii) reflect: Li2(1/(1-x)) = -Li2(x) - Li2(1-1/(1-x)) ... use: Li2(1/(1-x)) = z2? Use identity:
#     Li2(z) + Li2(1/(1-z))... simpler: Li2(1/(1-x)) = z2 - log(1-x)*log(x/(x-1)^2)? -> use Euler: Li2(w)+Li2(1-w)=z2-log(w)log(1-w), w=1/(1-x), 1-w=-x/(1-x)
#     log(w) = -log(1-x), log(1-w) = log(-x) - log(1-x)
# (iii) Li2(1+x): inversion-free; keep. Hmm R1 has Li2(1+x) (the d-rep). Function identity for Li2(1+x)? reflection: Li2(1+x)+Li2(-x)=z2-log(1+x)log(-x). But Li2(-x) isn't in the R1 expression... derivative may bring it.
sub_id = {
  sp.polylog(2, x2f): -sp.polylog(2, x) - sp.log(1-x)**2/2,
  sp.polylog(2, r2hf): sp.pi**2/6 - (-sp.log(1-x))*(sp.log(-x)-sp.log(1-x)) - sp.polylog(2, -x/(1-x)),
}
F2 = F1e.subs(sub_id, simultaneous=True)
F2 = sp.expand(F2)
print("after dilog function-identities, terms:", len(sp.Add.make_args(F2)))
for term in sp.Add.make_args(F2): print("   ", term)
