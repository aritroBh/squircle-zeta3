import subprocess
PY = '<project>/.venv/bin/python'
test = r'''
import pickle, sympy as sp
import mpmath as mp
mp.mp.dps = 60
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
K1F = '<project>/research/08_theory/K1f_hexagon_order5'
bb = pickle.load(open(f'{K1G}/build_blocks.pkl','rb'))
JLOC_n = bb['JLOC_n']
E5 = pickle.load(open(f'{K1F}/E5_polys.pkl','rb'))[5]
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
z2s, z3s, z5s = sp.symbols('z2 z3 z5')
target = mp.mpf("656.4153619676247544766476")
def evalP(lA0v, lA1v, lBv):
    rd = {A0: mp.mpf('0.5'), A1: mp.mpf('0.5'), B: mp.mpf('0.5'), kap: mp.mpf(0),
          L2: mp.mpf(1), q0: mp.mpf('-0.5'), V02: mp.mpf(1), lA0: lA0v, lA1: lA1v, lB: lBv}
    Earc = mp.mpf(0)
    for (n,a,b,c), co in E5.items():
        if c != 0: continue
        cov = co.subs(rd)
        if cov == 0: continue
        jv = JLOC_n[(n,a,b)].subs({z2s: mp.zeta(2), z3s: mp.zeta(3), z5s: mp.zeta(5)})
        Earc += mp.mpf(sp.N(cov, 40))*mp.mpf(sp.N(jv, 40))
    return 120*Earc*6
# convention candidates
l2, l3 = mp.log(2), mp.log(3)
cands = {
  'log(A/w)=log(sqrt3/2)': (mp.log(mp.sqrt(3)/2),)*3,
  'log(A)=log(1/2)': (mp.log(mp.mpf('0.5')),)*3,
  'log(A/w) w=1/2': (mp.log(mp.mpf(1)),)*3,
  'log(2A)': (mp.log(1)*0+l2,)*3,
  'log(sqrt3)': (mp.log(mp.sqrt(3)),)*3,
  'zero logs': (mp.mpf(0),)*3,
}
for name, logs in cands.items():
    v = evalP(*logs)
    print("%-28s P5 = %s   diff %s" % (name, mp.nstr(v, 16), mp.nstr(v-target, 4)))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=1800)
print(r.stdout); print(r.stderr[-600:])
