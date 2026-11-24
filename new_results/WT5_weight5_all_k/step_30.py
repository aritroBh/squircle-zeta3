import subprocess, hashlib
base='<project>'
def h(p):
    return hashlib.sha256(open(p,'rb').read()).hexdigest()[:12]
pairs=[('A3_W4_proof/w4_symbol.py','M5b_w6_closed_form/w4_symbol.py'),
       ('WT2_weight2_all_k/fastred.py','WT4_weight4_all_k/fastred.py'),
       ('WT2_weight2_all_k/binteg_10.pkl','WT4_weight4_all_k/binteg_10.pkl')]
for a,b in pairs:
    print(a, h(base+'/research/08_theory/'+a), ' | ', b, h(base+'/research/08_theory/'+b))
# head of w4_coproduct.py
with open(base+'/research/08_theory/A3_W4_proof/w4_coproduct.py') as f:
    t=f.read()
print(t[:1500])
