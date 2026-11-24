import flint, subprocess
print("flint version:", flint.__version__)
print([m for m in dir(flint.fmpz_mat) if 'solve' in m.lower() or 'dixon' in m.lower() or 'rref' in m.lower() or 'det' in m.lower()])
print(subprocess.run(['vm_stat'], capture_output=True, text=True).stdout[:600])