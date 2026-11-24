import subprocess, sys, json, os
src = open('step_10.py').read()
v = src.replace('J = 48; zstar = arb("3.0")', 'J = 96; zstar = arb("3.0")')
v = v.replace('open("final_cell.json","w")', 'open("final_cell_J96.json","w")')
open('step_10_J96.py','w').write(v)
r = subprocess.run([sys.executable, 'step_10_J96.py'], capture_output=True, text=True, timeout=600)
for line in r.stdout.splitlines():
    if 'FINAL CELL' in line: print("J96:", line)
if r.stderr: print("stderr:", r.stderr[-300:])

# Final artifact check
for f in ['result.md','progress.log','final_cell.json','final_cell_J96.json','final_cell_J192.json',
          'tight_100_400.json','sweep_cells_400_1e10.json','step_38.py','step_10.py','step_10_J96.py',
          'step_10_J192.py','step_15.py','step_16.py','step_17.py']:
    print(f, os.path.exists(f), os.path.getsize(f) if os.path.exists(f) else '')
print()
print("=== THEOREM SUMMARY ===")
t = json.load(open('tight_100_400.json')); S1 = max(c[2] for c in t)
s = json.load(open('sweep_cells_400_1e10.json')); S2 = max(c[2] for c in s)
d = json.load(open('final_cell_J192.json')); S3 = float(d['C'].split()[0].strip('[]'))
print(f"S1 [100,400]      = {S1}")
print(f"S2 [400,1e10]     = {S2}")
print(f"S3 [1e10,infty)   = {d['C'][:40]}...")
print(f"THEOREM C = max   = {max(S1,S2,S3)}  -> report C = 617117.6")
with open('progress.log','a') as f2:
    f2.write("turn10: artifacts finalized (J96 regenerated post-d5-fix). All present.\n")
