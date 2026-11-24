
import json
R = json.load(open('results_main.json'))
print("residuals (full) and scaled residuals:")
for k in sorted(R.keys()):
    num, pre, resid, scaled = R[k]
    print(k, '| resid =', resid, '| scaled =', scaled)
