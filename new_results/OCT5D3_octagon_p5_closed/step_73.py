txt = open("<project>/paper/ieee/main.tex").read()
i = txt.find("\\mathcal R&=(1-573")
j = txt.find("\\tag{T3}", i)
print(txt[i:i+2200])
