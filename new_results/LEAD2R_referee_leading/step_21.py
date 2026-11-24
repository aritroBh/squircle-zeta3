tex = open('../../../paper/squircle_note/main.tex').read()
i = tex.find('label{thm:lead}')
j = tex.rfind('begin{theorem', 0, i)
print(tex[j:i+5200])
