
import os, json
ROOT = '<project>'
# Check LEAD1 leading_coeffs.json and cw_dps420.json
for f in ['LEAD1_leading_genfun/leading_coeffs.json']:
    p = os.path.join(ROOT, 'research/08_theory', f)
    d = json.load(open(p))
    print(f, "keys:", list(d.keys())[:20])
    s = json.dumps(d)
    i = s.find('335687558224020357')
    print("  ...", s[max(0,i-80):i+80] if i>=0 else "not found")
