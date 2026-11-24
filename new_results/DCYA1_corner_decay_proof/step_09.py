import subprocess
# Read lines 330-455
print(subprocess.run(['sed', '-n', '330,455p', '<project>/research/08_theory/CX6_decay_rest/proof.md'], capture_output=True, text=True).stdout)
