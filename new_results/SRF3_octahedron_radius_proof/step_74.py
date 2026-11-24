
import os
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
r = open(wd+"/result.md").read()
print("has V-split:", "V-split" in r, "| has -sqrt(3pi/4) d-law:", "-\\sqrt{\\frac{3\\pi}4}\\,k^{-3/2}\\,\\mathrm{Im}[z_0" in r)
print("has pi-th/2=+2.2476162:", "+2.2476162" in r)
# print repaired proof sections for a final read
p = open(wd+"/proof.md").read()
i0 = p.find("### 3.2")
i1 = p.find("### 3.4")
print(p[i0:i1])
