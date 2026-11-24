
base='<project>/research/'
print(open(base+'08_theory/COR4_M3alpha_classical/step13_library.py').read()[:3000])
import json
sa=json.load(open(base+'08_theory/COR4_M3alpha_classical/sunit_args.json'))
print(type(sa), len(sa))
print(sa[:10] if isinstance(sa,list) else list(sa.items())[:10])
