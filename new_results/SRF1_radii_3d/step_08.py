
tex = open('<project>/paper/ieee/main.tex').read()
print(tex[75835:76200])
print("="*80)
# find lemma holo statement
i = tex.find('lem:holo', 75000)
print(tex[76000:79200])
