import subprocess
out = subprocess.run(['sed', '-n', '905,960p', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout
print(out)
