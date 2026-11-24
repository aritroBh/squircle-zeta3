import subprocess, os
root = '<project>'
# search for G1's corner formulas B2/B3
r = subprocess.run(['grep','-rl','B2', root + '/research/10_frontier', root + '/research/08_theory/D10c_corner_general'], capture_output=True, text=True)
print(r.stdout[:2000])
r2 = subprocess.run(['grep','-rn', "B2'", root + '/research/10_frontier/W1_corner_family'], capture_output=True, text=True)
print(r2.stdout[:1500])
