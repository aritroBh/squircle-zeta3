import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
import sympy as sp
u = sp.symbols('u')
S2 = sp.sqrt(2); S = sp.sqrt(1+u**2)
a = 1+S2; b = S2-1
v = u+S
meas = (1/(v+a)-1/(v-a))*sp.diff(v,u)
m_odd = (1+S2*u/S)/(1-u**2)
for uu in [0.3, -0.5, 0.87]:
    print(uu, sp.N(meas.subs(u,uu),25), sp.N(m_odd.subs(u,uu),25), sp.N((meas-m_odd).subs(u,uu),10))
