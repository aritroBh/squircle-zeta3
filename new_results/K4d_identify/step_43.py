import os
WORK='research/08_theory/K4d_identify'
print("=== bgA.log (h192?) ===")
print(open(WORK+'/bgA.log').read())
print("modelA_h192.json exists:", os.path.exists(WORK+'/modelA_h192.json'))
