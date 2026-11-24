import os, shutil, time
base = '<project>/research/08_theory'
tow3 = base + '/TOW3_tower_depth6'
work = base + '/TOW4_tower_depth8'
os.makedirs(work, exist_ok=True)

# copy canonical libs + data caches from TOW3 (never write into other job folders; copy what we reuse)
for f in ['tow2lib.py','tow3lib.py','L_parsed.pkl','L_simplified.pkl','arho56_simplified.pkl',
          'e61_closed.pkl','e_fits.pkl','rho_closed_forms.pkl','arho_simplified.pkl',
          'prove_V56.py','prove_adj56.py','prove_E3b56.py','prove_tops.py']:
    src = os.path.join(tow3, f)
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(work, f))
        print("copied", f, os.path.getsize(src))
    else:
        print("MISSING", f)

# also WTS kern files for Task C
wts = base + '/WTS_weight_k_structure'
print("WTS:", sorted(os.listdir(wts))[:30] if os.path.isdir(wts) else "no dir")
