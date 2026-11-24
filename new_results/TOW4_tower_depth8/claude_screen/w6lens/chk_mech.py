import sympy as sp
v=sp.symbols('v')
t=(v**2+2*v-1)/(4*v)
DL1=sp.diff(sp.log(t),v); DL2=sp.diff(sp.log(1-t),v)
print('t(1/v)-(1-t(v)) =', sp.simplify(t.subs(v,1/v)-(1-t)))
print('F(DL1)+DL2 =', sp.simplify(DL1.subs(v,1/v)/v**2 + DL2))
print('DL2 at v=i:', sp.simplify(DL2.subs(v,sp.I)), ' derivative at i nonzero:', sp.simplify(sp.diff(DL2,v).subs(v,sp.I)))
# compare with the reducer's tabulated DL: poles
print(sp.apart(sp.simplify(DL2),v))
