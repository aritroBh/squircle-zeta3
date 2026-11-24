import shutil, os
M4E='<project>/research/08_theory/M4e_w5_lchi'
WORK='<project>/research/08_theory/M4f_w5_extras'
files = ['w5_symbol.py','w5_search.py','w5_cobracket.py','w5_precision.py','nice_search.py',
         'nice_sat_s0.py','lchi_second.py','lchi_final.py','lchi_valrel.py','verify_nice.py',
         'display_nice.py','w5_verify.py',
         'W5_5400.json','Lchi8_5.json','cobracket.json','nice_ksat_s0.json','nice_symbol.json',
         'nice_atoms.pkl','lchi_final_1200.json','lchi_final_1500.json','closed_form_nice.md']
for f in files:
    shutil.copy(f'{M4E}/{f}', f'{WORK}/{f}')
os.makedirs(f'{WORK}/planner', exist_ok=True)
for f in ['R5_nice_reduced.json','planner_reduce.py','planner_eval_R5nice.py']:
    shutil.copy(f'{M4E}/planner/{f}', f'{WORK}/planner/{f}')
with open(f'{WORK}/progress.log','a') as fh:
    fh.write('step_03: copied M4e scripts+data into M4f (never editing M4e).\n')
print(sorted(os.listdir(WORK)))
