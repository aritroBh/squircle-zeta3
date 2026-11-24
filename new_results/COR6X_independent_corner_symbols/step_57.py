t = r'''
d = '<project>/research/10_frontier/G1_zeta3_any_angle/'
for f in ['result.md', 'verify.md', 'm3_hpl_check.py']:
    print("=====", f)
    try:
        print(open(d + f).read()[:4000])
    except Exception as e:
        print(e)
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout[:9000])
