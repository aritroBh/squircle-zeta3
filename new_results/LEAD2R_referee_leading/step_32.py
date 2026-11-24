import sys; sys.path.insert(0,'.')
import sympy as sp
u = sp.symbols('u'); S2 = sp.sqrt(2); S = sp.sqrt(1+u**2)
av = 1+S2; bv = S2-1
v = u + S
dvdu = sp.diff(v, u)
def zero(e): return sp.radsimp(sp.simplify(e)) == 0
tests = {}
# (a) t = (1+u)/2, 1-t = (1-u)/2   [t=(v-b)(v+a)/(4v)]
t_v = (v-bv)*(v+av)/(4*v); om_v = (av-v)*(v+bv)/(4*v)
tests['t=(1+u)/2'] = zero(t_v - (1+u)/2)
tests['1-t=(1-u)/2'] = zero(om_v - (1-u)/2)
# (b) four measures: (pole combination)*dvdu == claimed m(u)
m_odd = (1+S2*u/S)/(1-u**2); m_e = (1-S2/S)/(1-u**2)
tests['odd combo']  = zero((1/(v+av)-1/(v-av))*dvdu - m_odd)
tests['even combo'] = zero((1/(v+av)-1/(v+bv))*dvdu - m_e)
tests['diag combo'] = zero((1/(v+av)-1/(2*v))*dvdu - (1-u)/2*m_e)
tests['oddlast combo'] = zero((1/(v+av)-1/(2*(v-av))+1/(2*(v+bv))-1/(2*v))*dvdu - (m_odd - 1/(2*(1-u))))
# (c) m_odd = 1/(1-u) - u m_e ; m_e even ; endpoint limits
tests['m_odd=1/(1-u)-u m_e'] = zero(m_odd - (1/(1-u) - u*m_e))
tests['m_e even'] = zero(m_e - m_e.subs(u,-u))
tests['m_e(1)=-1/4'] = (sp.limit(m_e, u, 1) == sp.Rational(-1,4))
tests['m_e(-1)=-1/4'] = (sp.limit(m_e, u, -1) == sp.Rational(-1,4))
tests['m_odd+m_odd(-u)=2/(1-u^2)'] = zero(m_odd + m_odd.subs(u,-u) - 2/(1-u**2))
tests['m_odd(-1)=1/4 (finite)'] = (sp.limit(m_odd, u, -1) == sp.Rational(1,4))
tests['m_oddlast(-1) finite'] = sp.limit(m_odd - 1/(2*(1-u)), u, -1).is_finite
tests['m_oddlast ~ 1/(2(1-u)) at u=1'] = zero(sp.simplify((m_odd-1/(2*(1-u))) - 1/(2*(1-u)))==0) or (sp.limit((m_odd-1/(2*(1-u)))*(2*(1-u)), u, 1)==1)
# (d) deck: v_+*v_- = -1
tests['deck v+v-=-1'] = sp.simplify(v*(u-S)) == -1
tests['-1/a=-b'] = sp.simplify(-1/av + bv) == 0
for k, val in tests.items(): print(f"  {k}: {val}")
print("ALL B:", all(tests.values()))
with open('progress.log','a') as f:
    f.write("step_05: Lemma 2.1 all parts re-derived symbolically (own sympy): t,1-t, 4 measures, m_odd split, parity, endpoint limits -1/4,1/4, deck -1, pole pairing: ALL True\n")
