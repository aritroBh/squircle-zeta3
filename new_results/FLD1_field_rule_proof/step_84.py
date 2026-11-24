WD='<project>/research/08_theory/FLD1_field_rule_proof'
r=open(WD+'/result.md').read()
idx = r.find('k\\le5')
print(repr(r[idx-200:idx+80]))