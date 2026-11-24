import os
base = '<project>'
for p in ['research/08_theory', 'research/08_theory/FLD1_field_rule_proof']:
    full = os.path.join(base, p)
    print('===', p)
    for f in sorted(os.listdir(full)):
        st = os.stat(os.path.join(full, f))
        print(f'{st.st_size:9d}  {f}')
print('=== step_01.py (my cwd)')
print(open('step_01.py').read()[:2000])
