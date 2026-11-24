import os, glob
base = "<project>/research/08_theory/"
# find A3_W4_proof
for root, dirs, files in os.walk(base+"M4c_w5_closed_form"):
    for f in sorted(files):
        print(os.path.join(root,f).replace(base,""))
print("---- looking for A3_W4_proof ----")
cand = glob.glob("<project>/research/**/A3_W4_proof*", recursive=True)
print(cand[:20])
