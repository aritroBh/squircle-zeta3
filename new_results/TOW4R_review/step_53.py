
import os, re, glob
base = "<project>"
# paper
for p in glob.glob(base + "/paper/squircle_note/*.tex"):
    txt = open(p).read()
    if "w6k" in txt:
        for m in re.finditer(r"prop:w6k", txt):
            s = max(0, m.start()-1500)
            print(p)
            print(txt[s:m.start()+900])
            break
