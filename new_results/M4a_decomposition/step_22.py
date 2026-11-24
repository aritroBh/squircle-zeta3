import os, subprocess
BASE = '<project>'
TH = os.path.join(BASE,'research/08_theory')
# grep for numeric evaluation of words
for d in ['A3f3_level8_compute','A3h5_exact_certificate','A3h_weight4_target','A3f2_level8_w3w4','A3c_claude_elementary']:
    p = os.path.join(TH,d)
    if os.path.isdir(p):
        r = subprocess.run(['grep','-rln','polylog\|mpmath\|mp\.','--include=*.py',p],capture_output=True,text=True)
        print(d, ':', r.stdout[:800])
