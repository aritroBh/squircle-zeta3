import re
txt = open("<project>/paper/ieee/main.tex").read()
for pat in [r"C_m", r"local.{0,80}\\zeta", r"sin\^4", r"Z_e"]:
    ms = list(re.finditer(pat, txt))
    print(pat, len(ms))
# find the local-law definition
for m in re.finditer(r"C_m", txt):
    s = max(0, m.start()-900); print(txt[s:m.start()+1200]); print("="*90)
    break
