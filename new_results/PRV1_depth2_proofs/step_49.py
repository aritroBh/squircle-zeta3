wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()

old = """    checks = {"-x1/(1-x1) = x2": sp.simplify(-x1/(1-x1) - x2),
              "-k/(1-k) = -h": sp.simplify(-k/(1-k) + h),
              "log d = l/2 - L": sp.simplify(sp.log(d) - (l/2 - L)),
              "log k = -L": sp.simplify(sp.log(k) + L),
              "log x2 = -l/2 - L": sp.simplify(sp.log(x2) - (-l/2 - L)),
              "log(1-x1) = l/2": sp.simplify(sp.log(1-x1) - l/2)}
    for name, r in checks.items():
        print("     ", name, ": residual =", r)
        assert r == 0"""
new = """    # log identities are checked as exact MULTIPLICATIVE identities among the
    # positive real arguments (real log is injective on R_{>0}):
    checks = {"-x1/(1-x1) = x2": sp.simplify(-x1/(1-x1) - x2),
              "-k/(1-k) = -h": sp.simplify(-k/(1-k) + h),
              "d*(1+sqrt2) = sqrt2   [<= log d = l/2 - L]": sp.simplify(d*(1+s) - s),
              "k*(1+sqrt2) = 1       [<= log k = -L]": sp.simplify(k*(1+s) - 1),
              "x2*sqrt2*(1+sqrt2) = 1  [<= log x2 = -l/2 - L]": sp.simplify(x2*s*(1+s) - 1),
              "1-x1 = sqrt2          [<= log(1-x1) = l/2]": sp.simplify(1-x1 - s)}
    for name, r in checks.items():
        print("     ", name, ": residual =", r)
        assert r == 0"""
assert old in content
content = content.replace(old, new)
with open(f'{wd}/verify_full.py','w') as f: f.write(content)
print("part2b patched")

import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/_trial.py'], capture_output=True, text=True, cwd=wd, timeout=900)
print(r.stdout[-5000:])
print("STDERR:", r.stderr[-1500:])
