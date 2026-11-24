
import hashlib, os
def h(p): 
    with open(p,'rb') as f: return hashlib.md5(f.read()).hexdigest()
for f in ['poly_coeffs.py','step_09.py','step_10.py','z3coeff_m2_8.json','vals_m2_8_200.json']:
    mine = h(f) if os.path.exists(f) else 'MISSING'
    pred = h('../D15b_polygon_law/'+f) if os.path.exists('../D15b_polygon_law/'+f) else 'MISSING'
    print(f"{f:28s} mine={mine[:12]} pred={pred[:12]} {'SAME' if mine==pred else 'DIFF'}")
# check step_01 exists in mine
print("\nmy step_01.py exists:", os.path.exists('step_01.py'))
