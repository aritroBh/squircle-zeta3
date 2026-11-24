import subprocess
# Read full proof.md via sed in chunks - first get line count
out = subprocess.run(['wc', '-l', '<project>/research/08_theory/CX6_decay_rest/proof.md'], capture_output=True, text=True).stdout
print("Line count:", out)
# Read lines 90-210
print(subprocess.run(['sed', '-n', '90,200p', '<project>/research/08_theory/CX6_decay_rest/proof.md'], capture_output=True, text=True).stdout)
