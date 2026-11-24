
import os
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
log = open(wd+"/engine_K120.log").read()
print("LOG length:", len(log))
print(log[-500:])
print("files:", sorted(os.listdir(wd)))
