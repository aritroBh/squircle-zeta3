WORK='<project>/research/08_theory/M4f_w5_extras'
s = open(f'{WORK}/subset_exact.py').read()
s = s.replace("    rows = [M[i] + rhs[i] for i in range(86)]\n    mat, rank = flint.fmpq_mat(rows).rref()",
              "    rows = [[str(c) for c in M[i] + rhs[i]] for i in range(86)]\n    mat, rank = flint.fmpq_mat(rows).rref()")
s = s.replace("            Mx = [[coords[c][i] for c in allow] + [tbar[i]] for i in range(72)]\n            mat, rank = flint.fmpq_mat(Mx).rref()",
              "            Mx = [[str(coords[c][i]) for c in allow] + [str(tbar[i])] for i in range(72)]\n            mat, rank = flint.fmpq_mat(Mx).rref()")
open(f'{WORK}/subset_exact.py','w').write(s)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','subset_exact.py'],
                   cwd=WORK, capture_output=True, text=True, timeout=2400)
print('elapsed', time.time()-t0)
print(r.stdout)
print('STDERR', r.stderr[-2000:])
