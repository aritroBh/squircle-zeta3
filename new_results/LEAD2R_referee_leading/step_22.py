tex = open('../../../paper/squircle_note/main.tex').read()
i = tex.find('the $1/(1-u)$ part therefore contributes')
print(tex[i:i+3000])
