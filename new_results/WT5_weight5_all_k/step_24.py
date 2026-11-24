base='<project>'
with open(base+'/research/08_theory/M5b_w6_closed_form/w5_search.py') as fh:
    txt=fh.read()
i=txt.find('def numeric')
print(txt[i:i+4200])
