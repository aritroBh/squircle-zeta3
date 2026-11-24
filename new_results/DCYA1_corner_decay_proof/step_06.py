import subprocess
print("=== CX6_decay_rest/proof.md ===")
print(subprocess.run(['cat', '<project>/research/08_theory/CX6_decay_rest/proof.md'], capture_output=True, text=True).stdout)
