import os
base = "<project>/research/08_theory/"
D3 = base+"DOD3_dodecagon_symbol/"
for f in sorted(os.listdir(D3)):
    if f.startswith("step_") and f.endswith(".py"):
        txt = open(D3+f).read()
        if "R_li3_content" in txt or "R_all12_sym_str" in txt:
            print("=====", f, "len", len(txt))
