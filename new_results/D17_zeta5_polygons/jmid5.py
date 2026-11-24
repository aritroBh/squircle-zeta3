
import sympy as sp, sys, time, pickle
sys.path.insert(0, "research/08_theory/D17_zeta5_polygons")
import sym_D5
D, syms = sym_D5.build_D(5, reduced=True)
L0,L1 = syms['L0'], syms['L1']
s0,s1,c0,c1,S0,Sc = syms['s0'],syms['s1'],syms['c0'],syms['c1'],syms['S0'],syms['Sc']
D5 = sp.expand(D[5])
Gsum = sp.expand(D5.coeff(L0,3).coeff(L1,2) + D5.coeff(L0,2).coeff(L1,3))
sf,cf,Sh,Ch = sp.symbols('sf cf Sh Ch')
subsd = {sp.sqrt(1+Sc**2/S0**2): 1/cf, s0: sf*Ch+cf*Sh, c0: cf*Ch-sf*Sh,
  s1: Sh*cf-Ch*sf, c1: -(Ch*cf+Sh*sf), S0: cf/Sh, Sc: -sf/Sh}
g = sp.expand(Gsum.subs(subsd))
for pw in [2,3,4,5,6,7,8]:
    g = sp.expand(g.subs(sf**pw, sf**(pw-2)*(1-cf**2)))
terms = sp.Add.make_args(g)
even=[]
for t in terms:
    dd=t.as_powers_dict(); a=dd.get(sf,0); b=dd.get(cf,0); p=dd.get(Sh,0); q=dd.get(Ch,0)
    rest=sp.simplify(t/(sf**a*cf**b*Sh**p*Ch**q))
    if a==0: even.append((b,p,q,rest))
    else: print("ODD term!")
print("J_mid: %d even terms"%len(even))
def INT(b):
    if b==-2: return 2*sp.tan(h)
    if b==-4: return 2*(sp.tan(h)+sp.tan(h)**3/3)
    if b==-6: return 2*(sp.tan(h)+sp.Rational(2,3)*sp.tan(h)**3+sp.tan(h)**5/5)
h=sp.Symbol('h',positive=True)
I=0
for b,p,q,r in even: I+=r*Sh**p*Ch**q*INT(b)
I=sp.expand(I).subs(sp.tan(h),Sh/Ch); I=sp.expand(I)
Ir=sp.expand(I).subs(Ch**2,1-Sh**2); Ir=sp.expand(Ir)
for pw in [3,4,5,6,7,8]: Ir=sp.expand(Ir.subs(Ch**pw,Ch**(pw-2)*(1-Sh**2)))
Ir=sp.simplify(Ir)
print("J_mid(h) =", sp.factor(Ir))
# z2z3 coeff of P5: C_m^{(2,3)} = 240 m * (-12) * J_mid = -2880 m J_mid
print("C_m^{(2,3)} = -2880 m J_mid =", sp.factor(-2880*sp.Symbol('m')*Ir))
with open('Jmid_of_h.txt','w') as f: f.write(sp.srepr(sp.factor(Ir)))
