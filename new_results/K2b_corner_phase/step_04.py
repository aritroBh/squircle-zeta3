
with open(base+'mk_pi_2.tsv') as f:
    lines = f.readlines()
print(len(lines))
for ln in lines[:3] + lines[100:102] + lines[300:302]:
    print(repr(ln))
