base = 'research/08_theory/'
import subprocess
# grep for E3_general and y(t) parametrization across K1 steps
print(subprocess.run(['grep','-rl','E3_general\|V0\|y(t)\|edge', base+'K1_nonlocal_law/'], capture_output=True, text=True).stdout)
print("=== esum_engine.py ===")
with open(base+'K1_nonlocal_law/esum_engine.py') as f:
    print(f.read())
