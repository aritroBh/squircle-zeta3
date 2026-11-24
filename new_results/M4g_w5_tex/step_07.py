txt=open('../../../paper/ieee/main.tex').read()
print(len(txt))
import re
for kw in ['squircle','W_4','M_4','lambda','\\lam']:
    idxs=[m.start() for m in re.finditer(re.escape(kw),txt)]
    print(kw, len(idxs), idxs[:10])
