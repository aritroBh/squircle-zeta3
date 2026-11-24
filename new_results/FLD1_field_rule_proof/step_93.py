WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/proof.md').read()
i=s.find('## Lemma C'); j=s.find('## Lemma D')
print(s[i:j])