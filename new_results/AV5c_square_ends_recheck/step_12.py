tex = open('../../../paper/ieee/main.tex').read()
i = tex.find('tab:closed')
# find the table definition
j = tex.find('label{tab:closed}')
print(tex[j-6000:j+500])
