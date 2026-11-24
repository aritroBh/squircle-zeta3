import subprocess
for pat, name in [('mobius_from_triple','map enumeration'), ('permutations(PTS,4)','cross-ratio table'),
                  ('hpl_to_cmzv(w)','hpl conversion test'), ('li_atom(','li atom'),
                  ("'written hexc2_common'|'hexc2_common.py','w'",'common module'),
                  ('fast_cmzv','fast evaluator'), ('cert_w23.json','cert23'),
                  ('w4_y.pkl','w4 solve'), ('WEIGHT-4 EXACT CERTIFICATE','w4 verify'),
                  ('zeta(2) = pi^2/6','zeta2 check'), ('numeric_cmzv.py','mp evaluator')]:
    out = subprocess.run(['grep','-l',pat]+[f'step_{i:02d}.py' for i in range(1,75)],capture_output=True,text=True).stdout
    print(f'{name}: {out.strip().splitlines()[:4]}')
