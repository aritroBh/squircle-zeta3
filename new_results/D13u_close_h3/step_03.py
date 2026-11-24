p = '<project>/research/08_theory/D13t_review_fixes/result.md'
txt = open(p).read()
i = txt.find('## R5')
j = txt.find('## R6')
print(txt[i:j])
