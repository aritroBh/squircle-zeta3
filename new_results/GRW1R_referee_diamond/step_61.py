# Final integrity check: confirm result.md, verify_grw1r.py, progress.log exist and count step files
import os
steps = sorted(f for f in os.listdir('.') if f.startswith('step_') and f.endswith('.py'))
print(len(steps), "step files:", steps[:6], "...", steps[-3:])
for f in ['result.md','verify_grw1r.py','progress.log']:
    print(f, os.path.getsize(f), "bytes")
with open('progress.log','a') as f:
    f.write("done: FINAL report issued. Deliverables: result.md, verify_grw1r.py (52s), progress.log, step_01..36\n")
