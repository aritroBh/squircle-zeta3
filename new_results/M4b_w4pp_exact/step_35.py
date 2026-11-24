import glob, os
fs = sorted(glob.glob('<project>/research/08_theory/M4b_w4pp_exact/step_*.py'), key=os.path.getmtime)
for f in fs[-6:]:
    print(f, os.path.getsize(f))
print('--- newest content ---')
print(open(fs[-1]).read())
