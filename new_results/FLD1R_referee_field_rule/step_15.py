import re, os
tex = open('<project>/paper/ieee/main.tex').read()
print(len(tex))
for lab in ['conj:field','thm:flat','thm:M2alpha','prop:M3angle','thm:M4alpha']:
    idx = tex.find(lab)
    print('='*100)
    print('LABEL', lab, 'at', idx)
    # print surrounding env
    start = tex.rfind('\\begin', 0, idx)
    end = tex.find('\\end', idx)
    print(tex[start:end+60][:2200])
