import os
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
print(sorted(os.listdir(G1))[:40])
src = open(G1 + 'corner_series.py').read()
print(src[:3000])
