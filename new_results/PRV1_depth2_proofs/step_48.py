wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()

old = """    for name, e in ids.items():
        r = sp.factor(diff(e))
        print("   d/dt[", name, "] =", r)
        assert r == 0
    # DUP3's derivative uses DUP2 as a function identity: check that reduction too
    r = sp.factor(diff(Tt2 - 4*Tt - 4*Tm).subs(Dt2, 2*Dt+2*Dm))
    print("   d/dt[DUP3] after DUP2 =", r)
    assert r == 0"""
new = """    for name, e in ids.items():
        r = sp.factor(diff(e))
        print("   d/dt[", name, "] =", r)
        if not name.startswith("DUP3"):
            assert r == 0
    # DUP3's derivative uses DUP2 as a function identity: check that reduction too
    print("   [d/dt DUP3 = (2/t)*(DUP2), so d/dt DUP3 = 0 on functions; check:]")
    r = sp.factor(diff(Tt2 - 4*Tt - 4*Tm).subs(Dt2, 2*Dt+2*Dm))
    print("   d/dt[DUP3] after DUP2 =", r)
    assert r == 0"""
assert old in content
content = content.replace(old, new)
with open(f'{wd}/verify_full.py','w') as f: f.write(content)
print("part1 patched")

import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/_trial.py'], capture_output=True, text=True, cwd=wd, timeout=900)
print(r.stdout[-5500:])
print("STDERR:", r.stderr[-1500:])
