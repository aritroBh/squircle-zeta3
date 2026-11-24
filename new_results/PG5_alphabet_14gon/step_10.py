import subprocess
# find where the 12-gon coefficients come from in PGEN1
out = subprocess.run(["grep","-rn","coeffs12\|data12\|N3_12\|12gon\|dodec","../PGEN1_polygon_alphabet/","--include=step_*.py","-l"],capture_output=True,text=True).stdout
print(out)
