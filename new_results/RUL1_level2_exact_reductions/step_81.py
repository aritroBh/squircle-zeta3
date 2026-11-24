import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import pickle, sympy as sp
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
ge = pickle.load(open(f'{K1G}/golden_edges.pkl','rb'))
print(type(ge))
print(ge if not hasattr(ge, 'keys') else list(ge.keys()))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=300)
print(r.stdout[:3000]); print(r.stderr[-300:])
