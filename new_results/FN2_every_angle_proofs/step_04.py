import os
base = "<project>"
with open(os.path.join(base, "paper/ieee/main.tex")) as f:
    tex = f.read()
idx = tex.find("Ahmed's Open Problem 3")
print("found at", idx)
# find sectioning command before that
start = tex.rfind("\\section", 0, idx)
end = tex.find("\\section", idx)
print("section span:", start, end)
section = tex[start:end]
print(section[:6000])
