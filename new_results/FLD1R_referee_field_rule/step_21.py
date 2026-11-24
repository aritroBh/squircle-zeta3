tex = open('<project>/paper/ieee/main.tex').read()
i = tex.find('Proof idea', 100962)
print(tex[100962+1500:100962+1500+2600])
