import os
for root, dirs, files in os.walk("<project>/research/08_theory"):
    for f in files:
        if ("R2_rhombus" in root or "D11_rhombus" in root) and f.endswith(".md"):
            print(os.path.join(root, f))
