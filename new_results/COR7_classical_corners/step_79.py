import subprocess
out = subprocess.run(['grep', '-n', 'eq:R}', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout
print(out)
out = subprocess.run(['grep', '-n', 'tab:closed', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout
print(out)
