import os, subprocess
Z3 = 'Zeta' + '(3)'
ROOT = os.path.join('<home>/Downloads/Research', Z3)
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.chdir(WORK)
PY = os.path.join(ROOT, '.venv/bin/python')
src = open('step05_numeric_Rpp.py').read()
i = src.find("# readable print, Conjecture-4 style")
tail = '''# readable print, Conjecture-4 style (coefficient (c_rat + c_sq*sqrt2) printed exactly)
def coef_str(a, b):
    if not b: return f'{a}'
    if not a: return f'{b}*sqrt2'
    return f'({a}{"+" if b > 0 else "-"}{abs(b)}*sqrt2)'
mc.log("readable form (R'' = sum (c_rat + c_sq*sqrt2) * atom):")
lines = [f'   {"+" if (a > 0 or (a == 0 and b >= 0)) else "-"} {coef_str(abs(a), abs(b))}*{n}'
         if False else f'   ({coef_str(a, b)}) * {n}' for n, a, b in terms]
for chunk in [lines[i:i+2] for i in range(0, len(lines), 2)]:
    mc.log(' ; '.join(chunk))
mc.log(f'step05 done s={time.time()-t0:.1f}')
'''
src = src[:i] + tail
open('step05_numeric_Rpp.py','w').write(src)
r = subprocess.run([PY, 'step05_numeric_Rpp.py'], capture_output=True, text=True, cwd=WORK)
print('EXIT', r.returncode)
print('\n'.join(r.stdout.splitlines()[3:]))
print(r.stderr[-800:] if r.returncode else '')
