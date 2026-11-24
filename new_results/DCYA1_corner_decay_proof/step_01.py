import subprocess
# Explore the work directory and inputs
print("=== CWD ===")
print(subprocess.run(['pwd'], capture_output=True, text=True).stdout)
print("=== CWD contents ===")
print(subprocess.run(['ls', '-la'], capture_output=True, text=True).stdout)
print("=== research/08_theory contents ===")
print(subprocess.run(['ls', '<project>/research/08_theory/'], capture_output=True, text=True).stdout)
