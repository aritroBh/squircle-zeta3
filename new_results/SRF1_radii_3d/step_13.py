
import json, os
d16 = '<project>/research/08_theory/D16_pure_weight'
print("=== planner_run_n3.py ===")
print(open(os.path.join(d16,'planner_run_n3.py')).read()[:3000])
print("=== hi_tau_numeric_n2.json (keys, first entries) ===")
j = json.load(open(os.path.join(d16,'hi_tau_numeric_n2.json')))
print(type(j), list(j)[:5] if isinstance(j, dict) else len(j))
s = json.dumps(j)[:500]; print(s)
