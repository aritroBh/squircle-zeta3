p = "<project>/paper/ieee/main.tex"
txt = open(p).read()
i = txt.find(r"\label{prop:M3angle}")
j = txt.find(r"\end{proposition}", i)
print(txt[j:j+6000])
