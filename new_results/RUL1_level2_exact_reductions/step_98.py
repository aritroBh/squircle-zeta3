import pathlib, subprocess
PY = '<project>/.venv/bin/python'
src = pathlib.Path('verify_rul1.py').read_text()
src = src.replace("""z2s_, z3s_, z5s_ = sp.symbols('z2s z3s z5s')
Earc = sp.Integer(0)
for (n,a,b,c), co in E5.items():
    if c != 0: continue
    cov = co.subs(rd)
    if cov == 0: continue
    jv = JLOC_n[(n,a,b)].subs({z2s_: pi2/6, z3s_: z3, z5s_: z5})
    Earc += cov*jv""",
"""z2g, z3g, z5g = sp.symbols('z2 z3 z5')
Earc = sp.Integer(0)
for (n,a,b,c), co in E5.items():
    if c != 0: continue
    cov = co.subs(rd)
    if cov == 0: continue
    jv = JLOC_n[(n,a,b)].subs({z2g: pi2/6, z3g: z3, z5g: z5})
    Earc += cov*jv""")
pathlib.Path('verify_rul1.py').write_text(src)
# standalone test of the Target 3 part
t3 = r'''
import sys, time, pickle, json
sys.dont_write_bytecode = True
import mpmath as mp
import sympy as sp
T00 = time.time()
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
K1F = '<project>/research/08_theory/K1f_hexagon_order5'
bb = pickle.load(open(f'{K1G}/build_blocks.pkl','rb'))
JLOC_n = bb['JLOC_n']
E5 = pickle.load(open(f'{K1F}/E5_polys.pkl','rb'))[5]
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
l2, l3, pi2, z3, z5 = sp.symbols('l2 l3 pi2 z3 z5')
rd = {A0: sp.Rational(1,2), A1: sp.Rational(1,2), B: sp.Rational(1,2), kap: sp.Integer(0),
      L2: sp.Rational(1,3), q0: sp.Rational(-1,6), V02: sp.Rational(1,3), lA0: -l2, lA1: -l2, lB: -l2}
z2g, z3g, z5g = sp.symbols('z2 z3 z5')
Earc = sp.Integer(0)
for (n,a,b,c), co in E5.items():
    if c != 0: continue
    cov = co.subs(rd)
    if cov == 0: continue
    jv = JLOC_n[(n,a,b)].subs({z2g: pi2/6, z3g: z3, z5g: z5})
    Earc += cov*jv
P5 = sp.expand(720*Earc)
RHS = (30720*l2 - 4920*pi2 - 61440*l2**2 + 6240*pi2*l2 + 30720*l2**3 + 448*pi2**2
       + 5760*z3*l2 - 1920*pi2*l2**2 - 5120*l2**4 + 11340*z5 + 5580*pi2*z3/6
       - 123*pi2**2*l2 - 1440*z3*l2**2 + 160*pi2*l2**3 + 256*l2**5)
print("(128/sqrt3)P5 - RHS =", sp.expand(128/sp.sqrt(3)*P5 - RHS) == 0)
P5e = sp.expand(P5)
cz5r = sp.simplify(sum(t for t in sp.Add.make_args(P5e) if t.has(z5) and not t.has(z3) and not t.has(pi2) and not t.has(l2))/z5)
cpi2z3r = sp.simplify(sum(t for t in sp.Add.make_args(P5e) if t.has(z3) and t.has(pi2) and not t.has(z5) and not t.has(l2))/(pi2*z3))
print("z5 coord:", cz5r, " 6*(pi2 z3) coord:", sp.simplify(6*cpi2z3r))
mp.mp.dps = 130
import rul1num
NV2 = rul1num.named_values(130)
d17 = json.load(open('<project>/research/08_theory/D17_zeta5_polygons/m3_total.json'))['P5']
P5_num = sp.N(P5.subs({l2: NV2['log2'], pi2: NV2['pi2'], z3: NV2['z3'], z5: NV2['z5']}), 125)
dP = abs(mp.mpf(str(P5_num)) - mp.mpf(d17))
print("|P5 - D17| =", mp.nstr(dP, 3))
print("time %.0fs" % (time.time()-T00))
'''
r = subprocess.run(['nice',PY,'-c',t3], capture_output=True, text=True, timeout=1200)
print(r.stdout); print(r.stderr[-500:])
