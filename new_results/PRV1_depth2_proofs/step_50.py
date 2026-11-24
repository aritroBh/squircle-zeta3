wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()

old = """    # golden trilogs: 3-term Landen @ 1/phi^2, DUP3 @ 1/phi
    solt = sp.solve([L3a + L3b + L3c - (z3 - l*z2 + 5*l**3/6),
                     L3a - 4*L3b - 4*L3c], [L3a, L3b, L3c], dict=True)[0]
    print("   Lewin: Li3(1/phi^2) =", solt[L3a])
    print("   sum relation: Li3(1/phi)+Li3(-1/phi) =", sp.expand(solt[L3b] + solt[L3c]))
    S = sp.expand(solt[L3b] + solt[L3c])"""
new = """    # golden trilogs: 3-term Landen @ 1/phi^2, DUP3 @ 1/phi.
    # Two equations in L3a,L3b,L3c: they determine L3a (Lewin's value) and the
    # SUM L3b+L3c exactly; the difference L3b-L3c remains free (classical:
    # Li3(1/phi), Li3(-1/phi) are individually independent of z3, z2, l).
    solt = sp.solve([L3a + L3b + L3c - (z3 - l*z2 + 5*l**3/6),
                     L3a - 4*L3b - 4*L3c], [L3a, L3b], dict=True)[0]
    print("   Lewin: Li3(1/phi^2) =", sp.expand(solt[L3a]))
    S = sp.expand(solt[L3b] + L3c)
    print("   sum relation: Li3(1/phi)+Li3(-1/phi) =", S)"""
assert old in content
content = content.replace(old, new)
with open(f'{wd}/verify_full.py','w') as f: f.write(content)
print("part3c patched")

import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/_trial.py'], capture_output=True, text=True, cwd=wd, timeout=900)
print(r.stdout[-2500:])
print("STDERR:", r.stderr[-1500:])
