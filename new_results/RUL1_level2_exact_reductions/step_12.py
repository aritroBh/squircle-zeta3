import subprocess
print(subprocess.run(['bash','-c','sed -n "1270,1360p" "<project>/paper/ieee/main.tex"'], capture_output=True, text=True).stdout)
