import os, subprocess
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
PY = os.path.join(ROOT, '.venv/bin/python')
src = open('step05_numeric_Rpp.py').read()
i = src.find("# readable print, Conjecture-4 style")
tail = '''# readable print, Conjecture-4 style
def qs(c1, c2):
    """format c1 + c2*sqrt2 (both assumed as given signs)."""
    if not c2: return f'{c1}'
    if not c1: return f'({c2})*sqrt2' if abs(c2) != 1 else ('sqrt2' if c2 == 1 else '-sqrt2')
    return f'({c1}{"+" if c2 > 0 else "-"}{abs(c2)}*sqrt2)'
out = []
for n, a, b in terms:
    sgn = '-' if (a < 0 or (a == 0 and b < 0)) else '+'
    out.append((sgn, qs(abs(a), abs(b)), n))
mc.log("readable form (R'' = sum of signed terms):")
for chunk in [out[i:i+3] for i in range(0, len(out), 3)]:
    mc.log('   ' + ' '.join(f'{s} {c}*{n}' for s, c, n in chunk))
mc.log(f'step05 done s={time.time()-t0:.1f}')
'''
src = src[:i] + tail
open('step05_numeric_Rpp.py','w').write(src)
r = subprocess.run([PY, 'step05_numeric_Rpp.py'], capture_output=True, text=True, cwd=WORK)
print('EXIT', r.returncode)
print(r.stdout)
print(r.stderr[-1000:] if r.returncode else '')
