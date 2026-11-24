
import os
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
p = open(wd+"/proof.md").read()
# show the exact regions to repair
i1 = p.find("**The zeros of $Q$.**")
i2 = p.find("**Theorem B")
print(p[i1:i2])
print("="*80)
i3 = p.find("### 3.1")
i4 = p.find("### 3.4")
print(p[i3:i4][:5200])
