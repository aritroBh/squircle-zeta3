
tex=open('<project>/paper/ieee/main.tex').read()
i=tex.find('eq:R}')
j=tex.find('\\begin{equation}\\label{eq:R}')
if j<0:
    # search for label{eq:R}
    import re
    for m in re.finditer(r'label\{eq:R\}', tex):
        print(m.start())
else:
    print(j)
print(tex[j:j+2500])
