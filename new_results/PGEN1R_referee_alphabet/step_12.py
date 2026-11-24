import subprocess
# find where PGEN1 computed the decagon pair (a,b)
print(subprocess.run(["grep","-rln","875\|31.*sqrt.*320\|chi_5\|L(χ","--include=step_*.py","../PGEN1_polygon_alphabet"],capture_output=True,text=True).stdout)
