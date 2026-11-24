import os, time
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import sympy as sp
u, S = sp.symbols('u S')
S2 = sp.sqrt(2)
a = 1 + S2; b = S2 - 1
v = u + S   # positive branch of u = (v-1/v)/2
def red(expr):
    expr = sp.together(expr)
    num, den = sp.fraction(expr)
    rel = sp.Poly(S**2-(1+u**2), S, extension=S2)
    num = sp.Poly(sp.expand(num), S, extension=S2); den = sp.Poly(sp.expand(den), S, extension=S2)
    _, numr = sp.div(num, rel); _, denr = sp.div(den, rel)
    d0 = sp.expand(denr.as_expr()).coeff(S,0); d1 = sp.expand(denr.as_expr()).coeff(S,1)
    norm = sp.expand(d0**2 - d1**2*(1+u**2))
    res = sp.expand(numr.as_expr()*(d0-d1*S))
    return sp.simplify(res.coeff(S,0)/norm), sp.simplify(res.coeff(S,1)/norm)
dvdu = v/S
combos = {'odd': 1/(v+a)-1/(v-a), 'even': 1/(v+a)-1/(v+b),
          'diag': 1/(v+a)-1/(2*v), 'oddlast': 1/(v+a)-1/(2*(v-a))+1/(2*(v+b))-1/(2*v)}
# targets
m_odd  = (1+S2*u/S)/(1-u**2)
m_e    = (1-S2/S)/(1-u**2)
m_diag = (1-u)/2*m_e
m_ol   = m_odd - 1/(2*(1-u))
for name, Rv in combos.items():
    c0,c1 = red(Rv*dvdu)
    meas = c0 + c1*S
    tgt = {'odd':m_odd,'even':m_e,'diag':m_diag,'oddlast':m_ol}[name]
    diff = sp.simplify(sp.together(meas-tgt).subs(S, sp.sqrt(1+u**2)))
    print(f"{name}: measure == target ?  {sp.simplify(diff)==0}")
# also m_odd = 1/(1-u) - u m_e
chk = sp.simplify(sp.together(m_odd - (1/(1-u)-u*m_e)).subs(S, sp.sqrt(1+u**2)))
print("m_odd = 1/(1-u) - u m_e :", chk==0)
# t = (1+u)/2 check: t(v) = (v-b)(v+a)/(4v) == (1+u)/2
t_v = (v-b)*(v+a)/(4*v)
c0,c1 = red(t_v)
print("t(v) = (1+u)/2 :", sp.simplify(c0+c1*S-(1+u)/2)==0)
# 1-t = (1-u)/2
u1 = (a-v)*(v+b)/(4*v)
c0,c1 = red(u1)
print("1-t = (1-u)/2 :", sp.simplify(c0+c1*S-(1-u)/2)==0)
with open('progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} step: exact u-measures (odd/even/diag/oddlast) + t=(1+u)/2 via quotient ring — sympy\n")
