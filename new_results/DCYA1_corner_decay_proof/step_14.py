import subprocess
print("=== D13d_decay_theorem/result.md ===")
print(subprocess.run(['cat','<project>/research/08_theory/D13d_decay_theorem/result.md'], capture_output=True, text=True).stdout[:9000])
