
import os, signal, subprocess
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
# The copy I made is actually K=120 (engine had not rerun before crash). Rename properly.
os.rename(wd+"/SRF3_hi_coeffs_K60.json", wd+"/SRF3_hi_coeffs_K120.json")
# kill the redundant K=120 regeneration
try:
    os.kill(27130, signal.SIGTERM)
    print("killed redundant engine")
except Exception as e:
    print("kill:", e)
print(sorted(os.listdir(wd)))
