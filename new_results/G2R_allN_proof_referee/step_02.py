import os
base = "<project>/research/"
for p in ["10_frontier/W1_corner_family/B_taylor_pi4/b3_allN_proof.py",
          "08_theory/G2_M4alpha/planner_check_b4.py"]:
    fp = os.path.join(base, p)
    print(fp, os.path.exists(fp), os.path.getsize(fp) if os.path.exists(fp) else None)
# also check cwd existence
cwd = os.path.join(base, "08_theory/G2R_allN_proof_referee")
print("cwd:", cwd, os.path.exists(cwd))
if os.path.exists(cwd):
    print(os.listdir(cwd))
