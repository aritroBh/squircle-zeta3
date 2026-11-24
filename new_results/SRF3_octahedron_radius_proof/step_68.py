
import os
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
p = open(wd+"/proof.md").read()
i = p.find("**Theorem D.**")
print(p[i:i+2400])
print("="*90)
j = p.find("*Derivation")
print(p[j:j+900])
print("="*90)
k = p.find("### 4.4")
print(p[k:k+1400])
