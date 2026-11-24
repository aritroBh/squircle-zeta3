import shutil, os, subprocess, datetime
base = '<project>/research/08_theory/'
work = base+'TOW3R_review/'
d3 = base+'TOW3_tower_depth6/'
os.makedirs(work+'tow3_copy', exist_ok=True)
for f in ['tow3lib.py','tow2lib.py','prove_V56.py','prove_adj56.py','prove_E3b56.py',
          'L_parsed.pkl','L_simplified.pkl','arho_simplified.pkl','arho56_simplified.pkl',
          'e61_closed.pkl','e_fits.pkl','rho_closed_forms.pkl','e61_fit.pkl','e61_samples.pkl',
          'e61_fit_system.pkl','Sderiv_fits.pkl']:
    src = d3+f
    if os.path.exists(src):
        shutil.copy(src, work+'tow3_copy/'+f)
    else:
        print('MISSING', f)
print(sorted(os.listdir(work+'tow3_copy')))
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_01 found stale step_33 (mult bug, transcript 696 bad first run; fixed rerun uncommitted); copied libs to tow3_copy; ratio=a5 confirmed at (3,2) by my code\n")
PY='<project>/.venv/bin/python'
print(subprocess.run([PY,'-c','import sympy, mpmath; print(sympy.__version__, mpmath.__version__)'],capture_output=True,text=True).stdout)
