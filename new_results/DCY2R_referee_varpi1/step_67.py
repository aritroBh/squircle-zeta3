import os
open('progress.log','a').write("DONE: verdict PROVED AFTER REPAIRS; verify_dcy2r.py 132s PASS; result.md complete\n")
files = sorted(os.listdir('.'))
steps = [f for f in files if f.startswith('step_')]
print("cwd files:", [f for f in files if not f.startswith('step_')])
print(f"{len(steps)} step files: {steps[0]} .. {steps[-1]}")
print(open('progress.log').read())
