
tex=open('<project>/paper/ieee/main.tex').read()
i=tex.find('thm:W4R')
print(tex[i-100:i+5500])
