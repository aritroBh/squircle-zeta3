import re
tex = open("<project>/paper/ieee/main.tex").read()
print(len(tex))
for label in ["thm:P2", "thm:P2asym"]:
    i = tex.find("\\label{%s}" % label)
    if i == -1:
        print(label, "not found"); continue
    # find the begin{theorem} before
    j = tex.rfind("\\begin{theorem}", 0, i)
    k = tex.find("\\end{theorem}", i)
    print("="*30, label, "="*30)
    print(tex[j:k+15])
