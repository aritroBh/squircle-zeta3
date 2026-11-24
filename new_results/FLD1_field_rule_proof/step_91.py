WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/proof.md').read()
i=s.find('## Lemma A'); j=s.find('## Lemma B')
print(s[i:j])