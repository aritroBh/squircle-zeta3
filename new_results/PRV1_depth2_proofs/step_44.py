import subprocess, sys
PY = '<project>/.venv/bin/python'
wd = 'research/08_theory/PRV1_depth2_proofs'
runner = '''
import verify_full as V
V.part0()
V.part1()
V.part2()
V.part2_specialisation()
V.part3()
V.part3_landen3()
V.part3_specialisation()
print("TRIAL RUN PARTS 0-3c COMPLETE")
'''
with open(f'{wd}/_trial.py','w') as f: f.write(runner)
r = subprocess.run([PY, f'{wd}/_trial.py'], capture_output=True, text=True, cwd=wd, timeout=600)
print(r.stdout[-4000:])
print("STDERR:", r.stderr[-2000:])
