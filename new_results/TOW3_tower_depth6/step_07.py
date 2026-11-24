import shutil, os
base = '<project>/research/08_theory'
work = f'{base}/TOW3_tower_depth6'
# Copy reusable library + pickles from TOW2 (allowed: COPY what you reuse)
for f in ['tow2lib.py', 'arho_simplified.pkl', 'L_simplified.pkl', 'rho_closed_forms.pkl', 'e_fits.pkl', 'Sderiv_fits.pkl']:
    src = f'{base}/TOW2_tower_all_depths/{f}'
    if os.path.exists(src):
        shutil.copy(src, f'{work}/{f}')
        print('copied', f)
# also copy proof scripts for reference
for f in ['prove_S12.py','prove_S34.py','prove_tops.py','prove_adj34.py','prove_E3b34.py']:
    src = f'{base}/TOW2_tower_all_depths/{f}'
    if os.path.exists(src):
        shutil.copy(src, f'{work}/{f}')
print(sorted(os.listdir(work)))
