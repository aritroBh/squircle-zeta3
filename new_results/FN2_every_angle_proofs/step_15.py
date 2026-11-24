import os
base = "<project>"
with open(os.path.join(base, "paper/ieee/README_verification_map.md")) as f:
    rm = f.read()
# find the row about "every angle"
i = rm.find("every angle")
print(rm[max(0,i-2000):i+2500])
