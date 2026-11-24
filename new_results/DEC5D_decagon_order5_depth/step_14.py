
import re
tex = open("<project>/paper/ieee/main.tex").read()
# find prop:order5
i = tex.find("prop:order5")
print(tex[i-3000:i+3000])
