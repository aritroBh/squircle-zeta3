import os, subprocess
d17 = "research/08_theory/D17_zeta5_polygons"
venv_py = "./.venv/bin/python"

script = r'''
import mpmath as mp, json
mp.mp.dps = 110
print("dps:", mp.mp.dps)
d = json.load(open("m3_total.json"))
rem = mp.mpf(d["rem"])
print("rem =", mp.nstr(rem, 50))
sqrt3=mp.sqrt(3); pi=mp.pi; l2=mp.log(2); l3=mp.log(3); z3=mp.zeta(3); z2=mp.zeta(2); G=mp.catalan
# broad hexagon elementary basis (no z5, no z2z3)
basis = [mp.mpf(1), sqrt3, pi, pi*sqrt3, pi**2, pi**2*sqrt3, pi**3, pi**3*sqrt3,
         pi**4, pi**4*sqrt3, pi**5, pi**5*sqrt3,
         l2, l2*sqrt3, l3, l3*sqrt3, l2**2, l3**2, l2*l3,
         z3, z3*sqrt3, z3*pi, z3*pi**2, z2, z2*sqrt3, G, G*pi, G*sqrt3, G*pi**2]
bn = ['1','r3','pi','pir3','pi2','pi2r3','pi3','pi3r3','pi4','pi4r3','pi5','pi5r3',
      'l2','l2r3','l3','l3r3','l2^2','l3^2','l2l3',
      'z3','z3r3','z3pi','z3pi2','z2','z2r3','G','Gpi','Gr3','Gpi2']
v = [rem]+basis
rel = mp.pslq(mp.matrix(v), tol=mp.mpf('1e-95'), maxcoeff=10**12, maxsteps=20000)
print("PSLQ (broad):", rel)
if rel is not None and rel[0]!=0:
    rel=list(rel)
    terms = {bn[k]: -mp.mpf(rel[k+1])/rel[0] for k in range(len(basis)) if rel[k+1]!=0}
    print("rem =", {k: mp.nstr(val,10) for k,val in terms.items()})
# smaller: divide by sqrt3
print("\nrem/sqrt3 =", mp.nstr(rem/sqrt3, 40))
print("rem/3 =", mp.nstr(rem/3, 40))
'''.replace("D17DIR", d17)
with open(os.path.join(d17,"m3_pslq.py"),"w") as f: f.write(script)
r = subprocess.run([venv_py, os.path.join(d17,"m3_pslq.py")], capture_output=True, text=True, cwd=d17, timeout=600)
print(r.stdout); print("ERR:", r.stderr[-2000:])
