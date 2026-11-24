import os, subprocess
# Grep K1 folders for how r / kappa is computed from edge data
for d in ['K1_nonlocal_law','K1c_engine_predict','K1d_engine_fixed']:
    dd = 'research/08_theory/%s/'%d
    if not os.path.isdir(dd): continue
    # find files mentioning kappa=r-1 or A0*A1*B
    r = subprocess.run(['grep','-rln','kappa'], dd, capture_output=True, text=True)
    files = [f for f in r.stdout.split() if f.endswith('.py')][:8]
    print("###", d, "files with 'kappa':", [os.path.basename(f) for f in files])
