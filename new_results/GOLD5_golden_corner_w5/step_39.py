import subprocess
files = ['../G2_M4alpha/step_%02d.py' % i for i in range(1, 89)]
out = subprocess.run(['grep', '-ln', 'tanh\|quad', *files], capture_output=True, text=True)
print(out.stdout)
out = subprocess.run(['grep', '-n', 'def .*quad\|mp.quad\|tanh', '../G2_M4alpha/step_01.py', '../G2_M4alpha/step_02.py', '../G2_M4alpha/step_03.py', '../G2_M4alpha/step_04.py', '../G2_M4alpha/step_05.py'], capture_output=True, text=True)
print(out.stdout[:2000])
