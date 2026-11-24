
tex = open('<project>/paper/ieee/main.tex').read()
# Find the lemma holo statement itself
i = tex.find('begin{lemma}', 39000)
print(tex[39800:42000])
