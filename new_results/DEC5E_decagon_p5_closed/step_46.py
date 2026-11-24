tot = bs['tot']
# print a few entries
items = sorted(tot.items(), key=lambda kv: kv[0])
print("total words in S(W5/L):", len(items))
for w, c in items[:12]:
    print(w, c)
print('...')
sub = bs['sub_sym']
for k in sub:
    print(k, len(sub[k]))
