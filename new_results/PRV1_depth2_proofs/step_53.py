wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()

old = '''def dig(d):
    d = abs(d)
    return "exactly 0" if d == 0 else f"1e{mp.log10(d):.1f}".replace("e-0","e-").replace("e-","e-")'''
new = '''def dig(d):
    d = abs(d)
    if d == 0:
        return "exactly 0"
    e = int(mp.floor(mp.log10(d)))
    m = float(d / mp.power(10, e))
    return f"{m:.1f}e{e}"'''
assert old in content
content = content.replace(old, new)
with open(f'{wd}/verify_full.py','w') as f: f.write(content)
print("dig patched")

import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/verify_full.py'], capture_output=True, text=True, cwd=wd, timeout=1800)
out = r.stdout
# print only Part 4 onwards
i = out.find('PART 4')
print(out[i-100:])
print("STDERR:", r.stderr[-1500:])
with open(f'{wd}/verify_full.out','w') as f:
    f.write(out + "\nSTDERR:\n" + r.stderr)
