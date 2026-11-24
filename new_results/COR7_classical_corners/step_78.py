import subprocess
out = subprocess.run(['sed', '-n', '960,1000p', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout
print(out)
