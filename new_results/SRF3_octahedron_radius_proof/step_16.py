
import os
base = "<project>"
d = os.path.join(base,"research/08_theory/SRF1_radii_3d")
# look at fit steps 33-36
for f in ["step_33.py","step_34.py","step_35.py","step_36.py"]:
    src = open(os.path.join(d,f)).read()
    print("="*30, f, "="*30)
    print(src[:3000])
