import time, subprocess
time.sleep(90)
print(open('step09b.log').read())
print("---")
# check flint methods available
r = subprocess.run(['<project>/.venv/bin/python','-c',
 'import flint; print("fmpz nullspace:", hasattr(flint.fmpz_mat, "nullspace")); print("fmpq nullspace:", hasattr(flint.fmpq_mat, "nullspace")); print("fmpz rref:", hasattr(flint.fmpz_mat,"rref")); print("nmod nullspace:", hasattr(flint.nmod_mat,"nullspace"))'],
 capture_output=True, text=True, cwd='.')
print(r.stdout, r.stderr)
