import subprocess
print("=== K2_corner_decay/result.md (first 7000 chars) ===")
print(subprocess.run(['cat','<project>/research/08_theory/K2_corner_decay/result.md'], capture_output=True, text=True).stdout[:7000])
