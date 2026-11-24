
with open('<project>/paper/ieee/main.tex') as fh:
    tex = fh.read()
import re
i = tex.find('tab:verify')
print(tex[133007:136000])
