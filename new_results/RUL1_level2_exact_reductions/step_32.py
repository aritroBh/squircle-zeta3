import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1num
mp.mp.dps = 50
for M in [20, 40, 60, 100]:
    got = rul1num.zeta_stuple_num(((2,1),), dps=45, M=M)
    print("M=%3d z2 err %8s" % (M, mp.nstr(abs(got-mp.pi**2/6),3)))
for M in [20, 40, 60, 100]:
    got = rul1num.zeta_stuple_num(((1,-1),(1,-1)), dps=45, M=M)
    print("M=%3d z(1,1;-1,-1) err %8s" % (M, mp.nstr(abs(got-(mp.log(2)**2-mp.pi**2/6)/2),3)))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-500:])
