src = open('step_14.py').read()
i = src.find('cases = {')
j = src.find('for nm, (c, logs, li2args) in cases.items():')
newcases = '''cases = {}
H = mp.mpf(1)/2
def L(*args):
    return list(args)
cases['1/2'] = (mp.mpf(1)/2,
    [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log3^2', mp.log(3)**2), ('log2log3', mp.log(2)*mp.log(3))],
    [H, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, -H, -mp.mpf(1)/3, -mp.mpf(2)/3, -mp.mpf(1)/4, -mp.mpf(3)/4, mp.mpf(1)/9])
cases['3/5'] = (mp.mpf(3)/5,
    [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log5^2', mp.log(5)**2), ('log2log5', mp.log(2)*mp.log(5)), ('log3^2', mp.log(3)**2)],
    [H, mp.mpf(1)/4, mp.mpf(3)/4, mp.mpf(1)/5, mp.mpf(2)/5, mp.mpf(3)/5, mp.mpf(4)/5, -H, -mp.mpf(1)/4, -mp.mpf(1)/5, -mp.mpf(2)/5, -mp.mpf(3)/5, -mp.mpf(4)/5, mp.mpf(1)/16, mp.mpf(3)/16, mp.mpf(5)/16, mp.mpf(9)/16])
cases['1/3'] = (mp.mpf(1)/3,
    [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log3^2', mp.log(3)**2), ('log2log3', mp.log(2)*mp.log(3))],
    [H, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/4, mp.mpf(3)/4, -H, -mp.mpf(1)/3, -mp.mpf(2)/3, -mp.mpf(1)/4, -mp.mpf(3)/4, mp.mpf(1)/9, mp.mpf(1)/8, mp.mpf(3)/8, mp.mpf(1)/16])
cases['4/5'] = (mp.mpf(4)/5,
    [('pi^2', mp.pi**2), ('log2^2', mp.log(2)**2), ('log5^2', mp.log(5)**2), ('log2log5', mp.log(2)*mp.log(5)), ('log3^2', mp.log(3)**2)],
    [H, mp.mpf(1)/5, mp.mpf(2)/5, mp.mpf(3)/5, mp.mpf(4)/5, mp.mpf(1)/9, mp.mpf(4)/9, mp.mpf(5)/9, -H, -mp.mpf(1)/5, -mp.mpf(2)/5, -mp.mpf(3)/5, -mp.mpf(1)/9, -mp.mpf(4)/9, mp.mpf(1)/16, mp.mpf(1)/25, mp.mpf(1)/8, mp.mpf(3)/8])
'''
src = src[:i] + newcases + src[j:]
open('step_14.py','w').write(src)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_14.py'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
