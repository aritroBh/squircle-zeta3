import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c","""
import os
d = '<project>/research/10_frontier/G1_zeta3_any_angle'
print(os.listdir(d))
"""], capture_output=True, text=True)
print(r.stdout)
r2 = subprocess.run(["<project>/.venv/bin/python","-c","""
import os
d = '<project>/research/08_theory'
print(sorted(os.listdir(d)))
"""], capture_output=True, text=True)
print(r2.stdout)
