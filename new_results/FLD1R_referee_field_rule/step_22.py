tex = open('<project>/paper/ieee/main.tex').read()
i = tex.find('\\label{thm:flat}')
j = tex.find('Proof idea', i)
print(tex[j:j+1900])
