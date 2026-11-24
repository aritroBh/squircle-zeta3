import os

# First, explore the working directory and inputs
base = "<project>"
print("CWD:", os.getcwd())
print()
# List my work dir
workdir = "research/08_theory/FN2_every_angle_proofs"
os.makedirs(workdir, exist_ok=True)
print("Work dir contents:", os.listdir(workdir) if os.path.exists(workdir) else "creating")
print()
# List AH3 dir
ah3 = os.path.join(base, "research/08_theory/AH3_fn_any_angle")
print("AH3 dir:", os.listdir(ah3))
print()
p3 = os.path.join(base, "research/07_push/P3_fn_graph/claude")
print("P3 claude dir:", os.listdir(p3))
