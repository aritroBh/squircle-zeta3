
import json
R = json.load(open('results_main.json'))
# print full-precision 2-L and residuals for the report
for k in sorted(R.keys()):
    num, pre, resid, scaled = R[k]
    print(k, '|2-L =', num[:36], '|resid =', resid[:12], '|scaled =', scaled[:10])
