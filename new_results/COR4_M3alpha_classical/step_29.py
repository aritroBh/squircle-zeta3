import re
paper = "<project>/paper/squircle_note/"
tex = open(paper+"main.tex").read()
i = tex.find(r"\begin{theorem}\label{thm:W4}")
print(tex[i:i+1400])
