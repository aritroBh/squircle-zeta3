import shutil
work = '<project>/research/08_theory/M5a_next_coefficient'
shutil.copy(f'{work}/step_19.py', f'{work}/omega5_parse.py')
with open(f'{work}/progress.log', 'a') as f:
    f.write("step: results_summary.json written; all 6 tasks complete\n")
print(open(f'{work}/progress.log').read())