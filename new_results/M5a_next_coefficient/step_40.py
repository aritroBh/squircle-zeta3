work = '<project>/research/08_theory/M5a_next_coefficient'
src = open(f'{work}/w6_solve.py').read()
src = src.replace("M = flint.fmpq_mat([[atomvecs[a][r] for a in sub] for r in range(4)])",
                  "M = flint.fmpq_mat([[str(atomvecs[a][r]) for a in sub] for r in range(4)])")
src = src.replace("pr = Mi * flint.fmpq_mat([[t] for t in T_rat])",
                  "pr = Mi * flint.fmpq_mat([[str(t)] for t in T_rat])")
src = src.replace("qr = Mi * flint.fmpq_mat([[t] for t in T_sqr])",
                  "qr = Mi * flint.fmpq_mat([[str(t)] for t in T_sqr])")
open(f'{work}/w6_solve.py', 'w').write(src)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice', PY, '-B', 'w6_solve.py'], capture_output=True, text=True, cwd=work)
print("elapsed:", time.time()-t0)
print(r.stdout)
print("STDERR:", r.stderr[:2000])