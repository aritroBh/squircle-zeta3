import shutil, os
base = '<project>/research/08_theory/FLD1_field_rule_proof/'
for f in ['zsum.py','hpl.py','step_20_bijection.py','step_21_fastmu.py','step_22_musymbolic.py',
          'step_23_zsum.py','step_24_inspect.py','step_25_construct.py','mu_exact_K5_N100.json',
          'mu_exact_K6_N60.json','mu_zsum_K5.pkl','mu_zsum_K6.pkl','verify_fld1.py']:
    shutil.copy(base+f, f)
open('progress.log','a').write('setup: copied FLD1 sources (read-only) into referee dir\n')
# quick smoke: zsum self-tests (validates shift/stuffle/hpl_coeff against direct definitions)
import subprocess
r = subprocess.run(['nice','<project>/.venv/bin/python','zsum.py'],
                   capture_output=True, text=True, cwd='.')
print(r.stdout, r.stderr[-500:] if r.stderr else '')
open('progress.log','a').write('step02a: zsum.py self-tests (shift/stuffle/harm/hpl_coeff_at vs direct Z_val & coeffs_exact): '+r.stdout.strip()+'\n')
