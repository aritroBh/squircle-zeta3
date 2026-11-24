
with open('../D15b_polygon_law/step_10.py') as f:
    print(f.read())
print("="*60)
import json
with open('../D15b_polygon_law/z3coeff_m2_8.json') as f:
    z3 = json.load(f)
for k,v in z3.items():
    print(k, ":", v if isinstance(v,(int,float)) else (v[:80] if isinstance(v,str) else v))
