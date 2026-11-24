base = "<project>"
with open(base + "/paper/ieee/main.tex") as f:
    tex = f.read()
j = tex.find("lem:holo")
# find the lemma environment containing lem:holo
start = tex.rfind("\\begin{lemma}", 0, j)
end = tex.find("\\end{lemma}", j) + len("\\end{lemma}")
print(tex[start:end])
