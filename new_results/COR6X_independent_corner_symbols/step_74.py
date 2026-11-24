t = r'''
import json
tbl = json.load(open('<project>/research/08_theory/COR6_true_point_symbols/cor6_table.json'))
for e in tbl:
    if e['label'].startswith(('n=5','n=8','n=12','c=cos','n=4','n=3')):
        def nz(d):
            return sum(1 for v in d.values() if any(x != "0" for x in v))
        print(f"{e['label']:30s} true components: {nz(e['true_d22']):3d}   naive components: {nz(e['naive_d22']):3d}   keys total: {len(e['true_d22'])}")
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout)
