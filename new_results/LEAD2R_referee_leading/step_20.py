tex = open('../../../paper/squircle_note/main.tex').read()
i = tex.find('thm:lead')
print(tex[max(0,i-3500):i+3000])
