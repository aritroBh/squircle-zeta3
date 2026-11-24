import subprocess
print("=== step_01.py in CWD ===")
print(subprocess.run(['cat', 'step_01.py'], capture_output=True, text=True).stdout)
print("\n=== CX6_decay_rest contents ===")
print(subprocess.run(['ls', '-la', '<project>/research/08_theory/CX6_decay_rest/'], capture_output=True, text=True).stdout)
print("\n=== K2b_corner_phase contents ===")
print(subprocess.run(['ls', '-la', '<project>/research/08_theory/K2b_corner_phase/'], capture_output=True, text=True).stdout)
