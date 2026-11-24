import subprocess
print("=== K2b_corner_phase/result.md ===")
print(subprocess.run(['cat', '<project>/research/08_theory/K2b_corner_phase/result.md'], capture_output=True, text=True).stdout)
