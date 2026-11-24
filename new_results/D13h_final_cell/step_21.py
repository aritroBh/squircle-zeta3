import subprocess, sys
src = open('step_10.py').read()
for J, tag in [(192, 'J192')]:
    v = src.replace('J = 48; zstar = arb("3.0")', f'J = {J}; zstar = arb("3.0")')
    v = v.replace('open("final_cell.json","w")', f'open("final_cell_{tag}.json","w")')
    open(f'step_10_{tag}.py','w').write(v)
    r = subprocess.run([sys.executable, f'step_10_{tag}.py'], capture_output=True, text=True, timeout=600)
    for line in r.stdout.splitlines():
        if 'FINAL CELL' in line or 'S2z' in line or 'W2reg' in line or 'W3reg' in line:
            print(f"J={J}:", line)
    if r.stderr: print("stderr:", r.stderr[-300:])
with open('progress.log','a') as f:
    f.write("turn5: final cell J-refinement: J=48 C=57.0676, J=96 C=54.3068, J=192 see stdout. Junction box [1e10,1.05e10] step_38 C=43592 (wide-box artifact; final cell far tighter). Spot checks of sweep OK.\n")
