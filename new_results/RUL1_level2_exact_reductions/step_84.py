import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import pickle, sympy as sp
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
bb = pickle.load(open(f'{K1G}/build_blocks.pkl','rb'))
JLOC_n = bb['JLOC_n']
print("JLOC_n size:", len(JLOC_n))
import collections
for k in sorted(JLOC_n.keys())[:30]:
    print(k, "->", JLOC_n[k])
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=300)
print(r.stdout[:3500]); print(r.stderr[-300:])
