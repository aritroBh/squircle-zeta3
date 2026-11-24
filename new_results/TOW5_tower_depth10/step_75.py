import os, time
work = '<project>/research/08_theory/TOW5_tower_depth10'
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} VERIFY_TOW5: 11/11 lemmas PASS, 623s < 40min. Tasks A-D complete, E recorded.\n")
print(open(work+'/progress.log').read())
print("="*80)
# File inventory
for f in sorted(os.listdir(work)):
    if f.startswith('step_') or f=='progress.log': continue
    print(f, os.path.getsize(os.path.join(work,f)))
