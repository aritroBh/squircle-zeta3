import os
d = "../AH3_fn_any_angle"
print(os.listdir(d))
# look for a result/proof md
for f in os.listdir(d):
    if f.endswith(".md"):
        p = os.path.join(d,f)
        print("=====", f, os.path.getsize(p))
