t = r'''
import flint
print([m for m in dir(flint.fmpz_mat) if 'lll' in m.lower() or 'red' in m.lower()])
print([m for m in dir(flint) if 'lll' in m.lower()])
m = flint.fmpz_mat([[1,2,3],[4,5,6],[7,8,10]])
print([m for m in dir(flint.fmpz_mat) if not m.startswith('_')])
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
