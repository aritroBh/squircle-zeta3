import subprocess
out = subprocess.run(['grep', '-n', 'W4R', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout
print(out)
