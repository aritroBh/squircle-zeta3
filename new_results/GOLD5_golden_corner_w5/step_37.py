import subprocess, json
out = subprocess.run(['grep', '-ln', 'quad', '../G2_M4alpha/step_0*.py'], capture_output=True, text=True)
print(out.stdout)
d = json.load(open('../G2_M4alpha/flat_exact_N640_K4.json'))
print(type(d), len(d))
print(str(d)[:1500])
