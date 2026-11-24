import subprocess
root = '<project>'
r = subprocess.run(['find', root, '-name', 'mk_integrand*', '-o', '-name', 'ibp_reduce*'], capture_output=True, text=True)
print(r.stdout[:3000])
