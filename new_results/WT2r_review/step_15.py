import subprocess, os
root = '<project>'
r = subprocess.run(['grep', '-rln', 'def v_integrand\|def t_integrand\|mk_integrand', root+'/research'], capture_output=True, text=True)
print(r.stdout[:3000])
