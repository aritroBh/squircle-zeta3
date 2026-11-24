
import os
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
p = open(wd+"/proof.md").read()
i = p.find("**Theorem D.**")
print(p[i:i+1500])
print("="*100)
k = p.find("## Open items")
print(p[k:k+1800])
