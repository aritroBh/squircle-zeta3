base = "<project>"
# Find the sec:fn section in main.tex
with open(os.path.join(base, "paper/ieee/main.tex")) as f:
    tex = f.read()
print(len(tex))
# find the section
idx = tex.find("Graphs of $f^n")
print(idx)
idx2 = tex.find("sec:fn")
print(idx2)
