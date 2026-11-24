# verify_tow2r.py -- TOW2R consolidated referee check. Run in TOW2R_review/tow2_copy
# (copies of TOW2's scripts+libs) with project venv python. Prints PASS/FAIL per item.
import subprocess, sys, os, time
PY = sys.executable
def run(script, *args, timeout=1700):
    t0 = time.time()
    r = subprocess.run(['nice','-n','10',PY,script,*map(str,args)],
                       capture_output=True, text=True, timeout=timeout)
    ok = (r.returncode == 0)
    print(f'{"PASS" if ok else "FAIL"}  {script} {args}  ({time.time()-t0:.0f}s)')
    print('   ', r.stdout.strip().replace(chr(10),' | ')[-160:])
    return ok
ok = True
ok &= run('prove_S12.py')                 # S1,S4 lemmas, all (d,m)
ok &= run('prove_S34.py','3'); ok &= run('prove_S34.py','4')   # S''',S''''/S
ok &= run('prove_tops.py')                # top equations, all s
ok &= run('step_48.py')                   # V(3) non-adjacent
ok &= run('step_49.py')                   # V(4) non-adjacent (~6 min)
ok &= run('step_52.py')                   # slope (4,1) non-adjacent
for ds in [(3,0),(3,1),(4,0),(4,1),(4,2)]: ok &= run('prove_adj34.py',*ds)
ok &= run('prove_E3b34.py','3'); ok &= run('prove_E3b34.py','4')
print('ALL PASS' if ok else 'SOME FAIL')
# plus referee's own modules in this folder: tow2r_lib.py (first-principles Phi,
# Laurent extraction, T1); executed checks logged in progress.log:
#  S^l/S from genintegrand definition vs N_l/Pi_l: 172 cells d=4..12, 0 bad (l=1..4)
#  P_r Stirling tails: recurrence/anchors/values/divisibility, all True
#  Phi_p first principles vs binteg: key (4,0), delta 0..4, 5/5 exact
