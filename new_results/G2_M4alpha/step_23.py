base = "<project>/research/"
import os
# Check values files for M4
print(open(base + "10_frontier/G1_zeta3_any_angle/values_pi_2.txt").read()[:800])
print("="*80)
# M4e result.md first part
with open(base + "08_theory/M4e_w5_lchi/result.md") as f:
    t = f.read()
print(t[:3000])
