tex = open('../../../paper/ieee/main.tex').read()
j = tex.find('label{tab:closed}')
print(tex[j+2600:j+7000])
