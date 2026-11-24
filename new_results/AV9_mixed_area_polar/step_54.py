
import subprocess, os
# final progress.log entries
with open('progress.log','a') as f:
    f.write("step_14..17: c_9..c_12 exact and verified\n")
    f.write("DONE: all computations complete; report assembled\n")
print(open('progress.log').read())
print("=== files ===")
for f in sorted(os.listdir('.')):
    if f.endswith('.py') or f.endswith('.log'):
        print(f, os.path.getsize(f))
print("=== scratchpad ===")
for f in sorted(os.listdir('scratchpad')):
    print(f, os.path.getsize('scratchpad/'+f))
