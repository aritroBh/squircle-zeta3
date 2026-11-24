import pathlib, subprocess
PY = '<project>/.venv/bin/python'
src = pathlib.Path('rul1_step06.py').read_text()
old = src[src.index("# zeta(5), zeta(2)zeta(3) coefficients"):src.index("# numeric check")]
new = '''# zeta(5), zeta(2)zeta(3) coefficients
P5e = sp.expand(P5)
cz5 = sp.Integer(0)
for t in sp.Add.make_args(P5e):
    if t.has(z5) and not t.has(z3) and not t.has(pi2) and not t.has(l2) and not t.has(l3):
        cz5 += t
cz5 = sp.simplify(cz5/z5)
print("z5 coefficient of P5 =", cz5)
assert sp.simplify(cz5 - 2835*sp.sqrt(3)/32) == 0
cpi2z3 = sp.Integer(0)
for t in sp.Add.make_args(P5e):
    if t.has(z3) and t.has(pi2) and not t.has(z5) and not t.has(l2) and not t.has(l3):
        cpi2z3 += t
cpi2z3 = sp.simplify(cpi2z3/(pi2*z3))
print("pi2*z3 coefficient of P5 =", cpi2z3, "; as z2z3 = pi2 z3/6:", sp.simplify(6*cpi2z3))
assert sp.simplify(6*cpi2z3 - 1395*sp.sqrt(3)/32) == 0
log("TARGET 3b PASS: zeta5 = 2835 sqrt3/32, zeta2zeta3 = 1395 sqrt3/32 exactly")

'''
src = src[:src.index("# zeta(5), zeta(2)zeta(3) coefficients")] + new + src[src.index("# numeric check"):]
pathlib.Path('rul1_step06.py').write_text(src)
r = subprocess.run(['nice',PY,'rul1_step06.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout[-1400:]); print(r.stderr[-800:])
