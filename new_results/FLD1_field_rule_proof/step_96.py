WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/proof.md').read()
i=s.find('## Lemma ledger')
print(s[i:])