# verify_tow3r.py -- run in TOW3R_review/ with project venv python. PASS/FAIL per item.
import subprocess, sys, os, time
W = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable
def run(name, script, *args, cwd=None, timeout=1700):
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY,script,*map(str,args)], capture_output=True, text=True,
                       cwd=cwd or W+'/tow3_copy', timeout=timeout)
    ok = r.returncode==0
    print(f'{"PASS" if ok else "FAIL"}  {name}  ({time.time()-t0:.0f}s)')
    return ok
ok = True
ok &= run('L_parsed==L(r) symbolic r=1..6', W+'/step_04.py')          # ~10 s
ok &= run('rho(3..6) 2886 cells + Bell==displayed', W+'/step_05.py')  # ~20 s
ok &= run('tau table + conjugation rules', W+'/step_15.py')           # ~5 s
ok &= run('fits==lib; Newton coefs 1932+2476 cells', W+'/step_12.py') # ~15 s
ok &= run('diagonal ehat; empirical E3 (16 cases)', W+'/step_13.py')  # ~60 s
ok &= run('own T1 recursion 44443 fresh cells', W+'/step_14.py')      # ~30 s
for s in [(5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3)]:
    ok &= run(f'prove_V56 {s}', 'prove_V56.py', *s)                   # ~8 s total
ok &= run('prove_adj56 (5,2)', 'prove_adj56.py', 5, 2)                # fast representative
ok &= run('prove_adj56 (6,3)', 'prove_adj56.py', 6, 3)
ok &= run('prove_E3b56 5', 'prove_E3b56.py', 5)
ok &= run('prove_E3b56 6', 'prove_E3b56.py', 6)
ok &= run('prove_S12', 'prove_S12.py'); ok &= run('prove_S34 3','prove_S34.py',3)
ok &= run('prove_S34 4', 'prove_S34.py', 4); ok &= run('prove_tops','prove_tops.py')
# SLOW (documented, rerun green in this review; command lines):
#   prove_adj56.py 5 0 (169 s) | 5 1 (25 s) | 6 0 (570 s) | 6 1 (103 s) | 6 2 (14 s)
#   TOW2 base: step_48 V(3) 30 s, step_52 slope(4,1), step_49 V(4) 256 s, prove_adj34 all five
# NOT RUN BY REFEREE: Task C finite kernel checks (kern_5 vs cubic(6..9) 40/40; ipos/ineg k=6..9)
print('ALL PASS' if ok else 'SOME FAIL')
