import os, pickle, json
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
with open('progress.log','a') as f:
    f.write('step_07: setup; read theory/genfun.json, check_top.out, w4/w6 identify docs\n')
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
# inspect kern pkl structure
with open(os.path.join(theory, 'kern_10.pkl'),'rb') as f:
    K10 = pickle.load(f)
print(type(K10))
if isinstance(K10, dict):
    ks = list(K10.keys())
    print(len(ks), ks[:10])
