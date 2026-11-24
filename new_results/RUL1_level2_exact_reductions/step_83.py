import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import pickle, sympy as sp
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
K1F = '<project>/research/08_theory/K1f_hexagon_order5'
bb = pickle.load(open(f'{K1G}/build_blocks.pkl','rb'))
J5, JLOC_n = bb['J5'], bb['JLOC_n']
print("JLOC_n keys sample:", list(JLOC_n.items())[:6])
E5 = pickle.load(open(f'{K1F}/E5_polys.pkl','rb'))[5]
keys0 = [k for k in E5 if k[3] == 0]
print("c=0 keys:", len(keys0), "of", len(E5))
for k in keys0[:8]:
    print(k, "->", E5[k])
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=600)
print(r.stdout[:4500]); print(r.stderr[-400:])
