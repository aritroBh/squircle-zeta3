
base='<project>/research/'
import json
FS = json.load(open(base+'08_theory/COR4_M3alpha_classical/F_symbol.json'))
print(type(FS), len(FS))
if isinstance(FS, dict):
    ks = list(FS.keys())[:10]
    for k in ks: print(k, FS[k])
else:
    print(FS[:5])
print(open(base+'08_theory/COR4_M3alpha_classical/step11_F_symbol.py').read()[:4000])
