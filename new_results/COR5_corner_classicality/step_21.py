import os
paper = "<project>/paper/ieee/main.tex"
with open(paper) as f:
    txt = f.read()
i = txt.find("\\label{eq:R}")
j = txt.find("tab:closed")
print(txt[j-200:j+7000])
