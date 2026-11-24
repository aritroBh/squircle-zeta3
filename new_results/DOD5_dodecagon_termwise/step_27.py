import os, subprocess
base = "<project>/"
p = base+"paper/ieee/main.tex"
print(os.path.exists(p))
if os.path.exists(p):
    txt = open(p).read()
    import re
    for m in re.finditer(r"regular decagon", txt):
        s = max(0, m.start()-1500); print(txt[s:m.start()+2500]); print("="*80)
