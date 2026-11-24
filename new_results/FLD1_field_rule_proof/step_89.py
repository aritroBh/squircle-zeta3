WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/proof.md').read()
# print Lemma A and Lemma C and Lemma D sections with markers
import re
for tag in ['## Lemma A','## Lemma B','## Lemma C','## Lemma D']:
    i=s.find(tag)
    j=s.find('## ', i+5)
    print("#### SECTION", tag, "len", (j-i) if j>0 else len(s)-i)
print("-----LEMMA D-----")
i=s.find('## Lemma D'); j=s.find('## Lemma E')
print(s[i:j])