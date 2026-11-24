import subprocess
# Read lines 200-320
print(subprocess.run(['sed', '-n', '200,330p', '<project>/research/08_theory/CX6_decay_rest/proof.md'], capture_output=True, text=True).stdout)
