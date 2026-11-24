import subprocess
PY = '<project>/.venv/bin/python'
test = '''
import sys; sys.dont_write_bytecode=True
import mpmath as mp
import rul1num
mp.mp.dps = 50
import time
for dps in [100, 200, 260]:
    t=time.time()
    got = rul1num.zeta_stuple_num(((2,1),), dps=dps, M=400)
    print("dps=%3d z2 err %8s (%.1fs)" % (dps, mp.nstr(abs(got-mp.pi**2/6),3), time.time()-t))
for dps in [200, 260]:
    t=time.time()
    got = rul1num.zeta_stuple_num(((1,-1),(1,-1)), dps=dps, M=400)
    print("dps=%3d z(1,1;-1,-1) err %8s (%.1fs)" % (dps, mp.nstr(abs(got-(mp.log(2)**2-mp.pi**2/6)/2),3), time.time()-t))
    got = rul1num.zeta_stuple_num(((2,-1),(1,1)), dps=dps, M=400)
    print("dps=%3d z21bar err %8s (%.1fs)" % (dps, mp.nstr(abs(got-mp.zeta(3)/8),3), time.time()-t))
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=1200)
print(r.stdout); print(r.stderr[-500:])
