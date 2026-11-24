
import os, json
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
ROOT = '<project>'
# Get the 405-digit c(6) value
d = json.load(open(os.path.join(ROOT, 'research/08_theory/LEAD1_leading_genfun/leading_coeffs.json')))
c6_ref = d['values']['6']['value']
print("LEAD1 c(6) digits:", d['values']['6']['digits'])
print("value:", c6_ref[:80])
# save locally for verify script
json.dump({'digits': d['values']['6']['digits'], 'value': c6_ref}, open(os.path.join(WD, 'c6_lead_reference.json'), 'w'), indent=1)

import mpmath as mp
mp.mp.dps = 170
Evals = json.load(open(os.path.join(WD, 'e5m_values.json')))['values']
E55 = mp.mpf(Evals['E5,5'])
c6 = E55 / 24
diff = abs(c6 - mp.mpf(c6_ref))
print("E_5,5/24 - LEAD1 c(6) =", mp.nstr(diff, 3))
