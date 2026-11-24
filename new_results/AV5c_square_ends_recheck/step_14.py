tex = open('../../../paper/ieee/main.tex').read()
j = tex.find('label{tab:closed}')
k = tex.find('end{table*}', j)
seg = tex[j:k]
print(seg[:2900])
