import subprocess, time, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
PY = '<project>/.venv/bin/python'
code = '''
import sys, time
sys.path.insert(0, %r)
import sympy as sp
import tow5lib as T5, tow3lib as T3, tow2lib as T2
t0=time.time()
san1 = all(sp.cancel(sp.together(T5.forced(dl,s)-T3.e(dl,s)))==0 for dl in range(0,7) for s in range(0,dl//2+1))
print("F == tow3lib.e (delta<=6, all s):", san1, flush=True)
san2 = all(sp.cancel(sp.together(T5.forced(dl,dl//2)-T2.e_top(dl)))==0 for dl in range(0,11))
print("F(delta,top) == zigzag e_top (delta<=10):", san2, f"({round(time.time()-t0,1)}s)", flush=True)
# print the new closed forms F(9,s), F(10,s) tops for the record
for dl in [9,10]:
    print(f"F({dl},{dl//2}) [top] =", sp.cancel(T5.forced(dl, dl//2)), flush=True)
''' % (work,)
t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'-c',code], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print("ERR", r.stderr[-500:] if r.stderr else "", "elapsed", round(time.time()-t0,1))
