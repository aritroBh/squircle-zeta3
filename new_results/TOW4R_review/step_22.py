
import os, shutil, sys, time
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
tow4 = base + "/research/08_theory/TOW4_tower_depth8"
os.makedirs(work, exist_ok=True)
# scratch copy of TOW4 sources (read-only usage; never write into TOW4)
tc = os.path.join(work, "tow4_copy")
if os.path.exists(tc): shutil.rmtree(tc)
os.makedirs(tc)
for f in ["tow2lib.py","tow3lib.py","tow4lib.py","forced78.py","prove_adj78_fast.py","prove_adj78.py",
          "prove_E3b78.py","prove_E3b56.py","prove_adj56.py","check_full78.py","check_data.py",
          "verify_tow4.py","L_simplified.pkl","L_parsed.pkl","arho_simplified.pkl","arho56_simplified.pkl",
          "forced_F78.pkl","a7_closed.pkl","a8_closed.pkl","e_fits.pkl","e61_closed.pkl","prove_S34.py","prove_S12.py","prove_tops.py"]:
    src = os.path.join(tow4, f)
    if os.path.exists(src): shutil.copy(src, tc)
    else: print("MISSING", f)
# TOW2 originals for diffing
t2 = os.path.join(work, "tow2_orig"); os.makedirs(t2, exist_ok=True)
for f in ["tow2lib.py","prove_S34.py","prove_S12.py","prove_tops.py","L_simplified.pkl"]:
    src = base + "/research/08_theory/TOW2_tower_all_depths/" + f
    if os.path.exists(src): shutil.copy(src, t2)
# TOW3 originals
t3 = os.path.join(work, "tow3_orig"); os.makedirs(t3, exist_ok=True)
for f in ["tow3lib.py","prove_adj56.py","prove_E3b56.py","prove_V56.py","L_parsed.pkl"]:
    src = base + "/research/08_theory/TOW3_tower_depth6/" + f
    if os.path.exists(src): shutil.copy(src, t3)
    else: print("MISSING TOW3", f)
# progress log
with open(os.path.join(work, "progress.log"), "a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} setup: tow4_copy/tow2_orig/tow3_orig created; read result.md, proof_depth8.md, libs, CX9+TOW3R reports\n")
print(sorted(os.listdir(tc)))
