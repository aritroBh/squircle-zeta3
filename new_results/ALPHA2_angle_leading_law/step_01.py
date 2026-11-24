import subprocess
# Explore the directory structure
print(subprocess.run(['find', '<project>/research/08_theory', '-maxdepth', '1', '-type', 'd'], capture_output=True, text=True).stdout)
