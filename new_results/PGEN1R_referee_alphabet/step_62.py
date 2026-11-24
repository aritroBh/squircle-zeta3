import pickle, sympy as sp
d = pickle.load(open('../OCT6_octagon_closed/final_tables.pkl','rb'))
tot = d['tot']
li3 = {}
for w, nm, cf in tot:
    if w == 3 and nm.startswith('L3'):
        e = sp.expand(cf)
        li3[nm] = e
        print(nm, '=', e)
# PGEN1's printed delta lists
import subprocess, os
print(subprocess.run(["grep","-n","delta\\|∧","../PGEN1_polygon_alphabet/verify_output.txt"],capture_output=True,text=True).stdout[:200])
