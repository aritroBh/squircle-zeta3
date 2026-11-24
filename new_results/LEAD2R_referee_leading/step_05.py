txt = open('../LEAD2_leading_proof/proof_lead.md').read()
lines = txt.split('\n')
print(len(lines))
# print section headers with line numbers
for i, l in enumerate(lines):
    if l.startswith('#') or l.startswith('**L') or l.startswith('**T') or l.startswith('**M') or l.startswith('**C'):
        print(i, l[:110])
