WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/result.md').read()
i=s.find('## Labels')
j=s.find('## Checks')
print(s[i:j])