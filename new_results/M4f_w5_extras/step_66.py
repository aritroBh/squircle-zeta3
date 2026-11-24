WORK='<project>/research/08_theory/M4f_w5_extras'
txt = open(f'{WORK}/w5_closed_form.tex').read()
i = txt.find(r'\operatorname{Li}_5')
print(txt[i-60:i+1700])
