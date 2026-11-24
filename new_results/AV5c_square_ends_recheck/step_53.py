tex = open('../../../paper/ieee/main.tex').read()
i = tex.find('label{tab:corner}')
j = tex.find('end{table', i)
seg = tex[i:j]
print(seg[:5500])
