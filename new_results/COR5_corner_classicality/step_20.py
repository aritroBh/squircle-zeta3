import os, re
paper = "<project>/paper/ieee/main.tex"
with open(paper) as f:
    txt = f.read()
i = txt.find("thm:W4R")
print(txt[max(0,i-3500):i+4500])
