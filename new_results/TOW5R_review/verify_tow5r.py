"""verify_tow5r.py -- TOW5R referee runner (< 10 min core). Run with the project venv python
from research/08_theory/TOW5R_review/:   PY verify_tow5r.py
Every item PASS/FAIL. Heavy certificates (already run by the referee; commands printed at the
end): verify_L.py (214 s), step_05.py full (131 s), step_14.py TOW5-data 2x425 (64 s),
step_15.py fresh rho d=41,42 (64 s), step_09.py sanity+seed-blindness (75 s),
step_18.py fresh towers deg 19/20 (25 s), step_19b.py MASTER (72 s), taskE_defect.py (70 s),
verify_tow5.py full (618 s, rerun by referee: 11/11).
"""
import subprocess, sys, os, time
W = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(W, 'tow5_copy')
PY = sys.executable
RESULTS = []
def run(name, script, *args, cwd=None, timeout=1700, must=None, expect_fail=False):
    t0 = time.time()
    r = subprocess.run(['nice', '-n', '10', PY, script, *map(str, args)],
                       capture_output=True, text=True, cwd=cwd or W, timeout=timeout)
    ok = (r.returncode != 0) if expect_fail else (r.returncode == 0)
    if must:
        ok = ok and all(m in r.stdout for m in must)
    RESULTS.append((name, bool(ok)))
    print(f'[{"PASS" if ok else "FAIL"}] {name}  ({round(time.time()-t0,1)}s)', flush=True)
    if not ok:
        print('  STDOUT tail:', r.stdout[-300:]); print('  STDERR tail:', r.stderr[-300:])
    return ok

t00 = time.time()
# 1. tau(9),(10) exact + 40-digit
run('TAU910 exact QSI + mp', 'step_06.py')
# 2. R_5 proof rerun
run('R5 prove_S5 rerun', 'prove_S5.py', cwd=C, must=['PROVED for all (d,m)'])
# 3. R_5 independent Leibniz re-derivation (126/126 terms, cancel==0)
run('R5 independent Leibniz', 'step_04.py', must=['STEP_04 VERDICT'])
# 4. R_5 own S_m(t) direct 5th derivative, 45 keys
run('R5 own S direct (quick)', 'step_05q.py', must=['bad = 0'])
# 5. L(9),L(10) from scratch (own S, own Q(i) series), 110 values
run('L(9,10) from scratch', 'step_07b.py', must=['bad = 0'])
# 6. closing F(9,5)=0, F(10,6)=0 rebuilt from scratch (delete the cache first)
pkl = os.path.join(C, 'forced_F910.pkl')
if os.path.exists(pkl): os.remove(pkl)
run('CLOSE rebuild delta=9', 'forced910.py', '9', cwd=C, must=['CLOSING delta=9: F(9,5) numerator == 0 : True'])
run('CLOSE rebuild delta=10', 'forced910.py', '10', cwd=C, must=['CLOSING delta=10: F(10,6) numerator == 0 : True'])
# 7. all 11 adjacent identities
for dl, s in [(9,0),(9,1),(9,2),(9,3),(9,4),(10,0),(10,1),(10,2),(10,3),(10,4),(10,5)]:
    run(f'ADJ({dl},{s})', 'prove_adj910_fast.py', dl, s, cwd=C, must=['numerator == 0 : True'])
# 8. fast-vs-slow + mutation
run('ADJ fast-vs-slow', 'step_11.py', must=['STEP_11 VERDICT: PASS'])
src = open(os.path.join(C, 'prove_adj910_fast.py')).read()
mut = src.replace('(-i)**(-4-t) * binom * i**(-6-delta)', '(-i)**(-4-t) * binom * i**(-4-delta)')
open(os.path.join(C, 'prove_adj910_MUT.py'), 'w').write(mut)
run('ADJ mutation must FAIL', 'prove_adj910_MUT.py', 9, 2, cwd=C, expect_fail=True)
os.remove(os.path.join(C, 'prove_adj910_MUT.py'))
# 9. E3 blindness
run('E3(9)', 'prove_E3b910.py', 9, cwd=C, must=['numerator == 0 : True'])
run('E3(10)', 'prove_E3b910.py', 10, cwd=C, must=['numerator == 0 : True'])
# 10. rho(9),(10) vs WTS binteg ground truth
run('RHO WTS binteg', 'step_22.py', must=['STEP_22 VERDICT: PASS'])
# 11. LAW vs TOW5 data (990+11 and 16225+66 cells)
run('LAW check_full910', 'check_full910.py', cwd=C, must=['bad = 0'])
run('LAW check_full910_ext', 'check_full910_ext.py', cwd=C, must=['bad = 0'])
# 12. LAW independent (WTS tower ground truth, deg 7..12, k<=14)
run('LAW WTS ground truth', 'step_17.py', must=['bad = 0'])
# 13. W7 (kernels) + own + LEAD beta_6 table
run('W7 official', 'verify_W7.py', cwd=C, must=['VERDICT W7: PASS'])
run('W7 own (incl. fresh k=14)', 'step_13b.py', must=['STEP_13B VERDICT: PASS'])
run('LEAD beta_6 table', 'step_21.py', must=['bad = 0'])

print('=' * 70)
npass = sum(1 for _, o in RESULTS if o)
print(f'VERIFY_TOW5R: {npass}/{len(RESULTS)} items PASS, total {round(time.time()-t00,1)}s', flush=True)
for nm, o in RESULTS: print(f'  {nm}: {"PASS" if o else "FAIL"}')
print()
print('Slow certificates (referee already ran; rerun commands):')
print('  cd tow5_copy && PY verify_L.py            # L cache == tow2lib L(r), r=1..10, symbolic (214 s)')
print('  PY step_05.py                             # own S_m d^5 direct, 180 keys (131 s)')
print('  PY step_14.py                             # rho(9),(10) vs TOW5 data 2x425 + binteg (64 s)')
print('  PY step_15.py                             # fresh rho d=41,42 vs a_9,a_10 (64 s)')
print('  PY step_09.py                             # CLOSE sanity + seed blindness dl=10 (75 s)')
print('  PY step_18.py                             # fresh towers deg 19,20, 1254+11 cells (25 s)')
print('  PY step_19b.py                            # MASTER x=-1 at 9,10 (72 s)')
print('  cd tow5_copy && PY taskE_defect.py        # defect prefactor law dl=2..10 (70 s)')
print('  cd tow5_copy && PY verify_tow5.py         # official runner, 11/11 (618 s)')
