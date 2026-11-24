G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
import os
print([f for f in os.listdir(G1) if 'flat' in f or f.endswith('.py')][:40])
