import os
base = "<project>"
with open(os.path.join(base, "paper/ieee/main.tex")) as f:
    tex = f.read()
start = tex.rfind("\\section", 0, 138319)
end = tex.find("\\section", 138319)
section = tex[start:end]
print("span:", start, end, "len:", len(section))
with open("<project>/research/08_theory/FN2_every_angle_proofs/sec_fn_extract.tex","w") as f:
    f.write(section)
print(section[:8000])
