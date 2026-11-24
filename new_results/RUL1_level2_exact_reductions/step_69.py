import pathlib, subprocess
PY = '<project>/.venv/bin/python'
src = pathlib.Path('rul1_step04b.py').read_text()
old = src[src.index("# expected closed forms"):src.index("d4 = sp.expand")]
new = '''# expected closed forms (D12 result.md sec 1.3); pi^4 = pi2^2/90
pi4 = pi2**2/90
N4_exp = -l2**5/120 + pi2/24*(1-l2+l2**2/2) + z3/4*(1-l2) + pi4/2880
N5_exp = (l2**6/720 + pi2/24*(1-l2+l2**2/2-l2**3/6) + z3/4*(1-l2+l2**2/2)
          + pi4/2880*(1-l2) - 7*pi4/1920 - z3*pi2/48 + z5/8)
'''
src = src[:src.index("# expected closed forms")] + new + src[src.index("d4 = sp.expand"):]
pathlib.Path('rul1_step04b.py').write_text(src)
r = subprocess.run(['nice',PY,'rul1_step04b.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout[-1800:]); print(r.stderr[-1200:])
