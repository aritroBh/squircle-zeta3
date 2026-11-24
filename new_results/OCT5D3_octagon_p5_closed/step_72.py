txt = open("<project>/paper/ieee/main.tex").read()
i = txt.find("eq:W3p")
print(txt[i-2500:i+500])
