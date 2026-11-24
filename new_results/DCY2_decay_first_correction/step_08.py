base = "<project>/"
txt = open(base+"paper/ieee/main.tex").read()
i = txt.find("sec:twoends")
# find the section label definition
j = txt.find("label{sec:twoends}")
print(j)
print(txt[j-2000:j+6000])
