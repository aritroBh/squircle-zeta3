
simp4 = '''
import sympy as sp
import sym_G
G21, G12, syms = sym_G.build()
s0,s1,c0,c1,S0,Sc = syms['s0'],syms['s1'],syms['c0'],syms['c1'],syms['S0'],syms['Sc']

phi, h = sp.symbols('phi h', positive=True)
subs = {
  s0: sp.sin(phi+h), c0: sp.cos(phi+h),
  s1: sp.sin(h-phi), c1: -sp.cos(h-phi),
  S0: sp.cos(phi)/sp.sin(h), Sc: -sp.sin(phi)/sp.sin(h),
}
def do(G,label):
    Ge = sp.expand(G).subs(sp.sqrt(1+Sc**2/S0**2), sp.sqrt(S0**2+Sc**2)/S0)
    Ge = Ge.subs(subs)
    Ge = sp.simplify(Ge)
    print("====",label,"====")
    print(sp.factor(Ge))
    return Ge
g21 = do(G21,'G21(phi)')
g12 = do(G12,'G12(phi)')
# check g12 = g21(-phi)
g21m = sp.simplify(g21.subs(phi,-phi))
print("G12 - G21(-phi) =", sp.simplify(g12-g21m))
'''
with open('step_05_subst.py','w') as f: f.write(simp4)
import subprocess, time
t0=time.time()
r = subprocess.run(['./.venv/bin/python','step_05_subst.py'],capture_output=True,text=True,timeout=600)
print("time",time.time()-t0)
print(r.stdout); print(r.stderr[-2000:] if r.stderr else "")
