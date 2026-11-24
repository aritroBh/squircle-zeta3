import subprocess
print("=== CX6_decay_rest/common.py ===")
print(subprocess.run(['cat','<project>/research/08_theory/CX6_decay_rest/common.py'], capture_output=True, text=True).stdout)
print("\n=== CX6_decay_rest/s2_assemble.py ===")
print(subprocess.run(['cat','<project>/research/08_theory/CX6_decay_rest/s2_assemble.py'], capture_output=True, text=True).stdout)
