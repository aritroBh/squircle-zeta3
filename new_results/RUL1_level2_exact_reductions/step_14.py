import subprocess
print(subprocess.run(['bash','-c','sed -n "1605,1635p" "<project>/paper/ieee/main.tex"'], capture_output=True, text=True).stdout)
