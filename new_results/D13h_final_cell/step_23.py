import subprocess, sys, json
src = open('step_10.py').read()
for J, tag in [(192, 'J192')]:
    v = src.replace('J = 48; zstar = arb("3.0")', f'J = {J}; zstar = arb("3.0")')
    v = v.replace('open("final_cell.json","w")', f'open("final_cell_{tag}.json","w")')
    open(f'step_10_{tag}.py','w').write(v)
    r = subprocess.run([sys.executable, f'step_10_{tag}.py'], capture_output=True, text=True, timeout=600)
    for line in r.stdout.splitlines():
        if 'FINAL CELL' in line or 'W3reg' in line:
            print(f"J={J}:", line)
    if r.stderr: print("stderr:", r.stderr[-300:])
d = json.load(open('final_cell_J192.json'))
print("J192 C =", d['C'])
with open('progress.log','a') as f:
    f.write(f"turn6: tail-q fix (valid joint bound exp<=e^0.665); J=48 C=53.6157, J=192 C={d['C']}.\n")
