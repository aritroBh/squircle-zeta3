import subprocess
print("=== D13s_fix_recertify/result.md ===")
print(subprocess.run(['cat','<project>/research/08_theory/D13s_fix_recertify/result.md'], capture_output=True, text=True).stdout[:8500])
