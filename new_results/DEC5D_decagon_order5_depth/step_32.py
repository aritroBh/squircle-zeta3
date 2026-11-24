
tex = open("<project>/paper/ieee/main.tex").read()
i = tex.find("thm:Cm")
print(tex[i-1500:i+1500])
