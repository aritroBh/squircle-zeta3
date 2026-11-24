import shutil, os, json, subprocess, time
base = '<project>/research/08_theory'
cwd = os.getcwd()
print("CWD:", cwd)
# Copy libraries (COPY, never edit in place)
for src in [f'{base}/ZQ5b_zagier_w5_larger_S/w5blib.py',
            f'{base}/ZQ5b_zagier_w5_larger_S/w5pipe.py',
            f'{base}/ZQ5_zagier_weight5/planner/planner_li5.py']:
    shutil.copy(src, cwd)
    print("copied", os.path.basename(src))
# Copy data dirs
os.makedirs('out', exist_ok=True)
for dd in ['d29','d13','d21','d5']:
    srcd = f'{base}/ZQ5b_zagier_w5_larger_S/out/{dd}'
    dstd = f'{cwd}/out/{dd}'
    os.makedirs(dstd, exist_ok=True)
    for f in os.listdir(srcd):
        if f.endswith('.json'):
            shutil.copy(os.path.join(srcd,f), dstd)
    print("copied out/"+dd, os.listdir(dstd))
# also d17_final.json as format reference
shutil.copy(f'{base}/ZQ5b_zagier_w5_larger_S/d17_final.json', cwd)
# start progress.log
with open('progress.log','w') as f:
    f.write(f"step_00_setup: libs copied (w5blib,w5pipe,planner_li5); data d29/d13/d21/d5; d17_final ref; {time.strftime('%H:%M:%S')}\n")
print(open('progress.log').read())
