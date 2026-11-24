import re
txt = open("<project>/paper/ieee/main.tex").read()
# find prop:order5
i = txt.find("prop:order5")
print(txt[i-200:i+3500])
