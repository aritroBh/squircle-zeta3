import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import pickle, sympy as sp
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
ge = pickle.load(open(f'{K1G}/golden_edges.pkl','rb'))
print("V:", ge["V"])
print("w:", ge["w"])
print("k:", ge["k"])
for i, e in enumerate(ge["edges"]):
    print("edge", i, e)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=300)
print(r.stdout[:4000]); print(r.stderr[-300:])
