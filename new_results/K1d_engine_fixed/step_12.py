import pickle
E3 = pickle.load(open('E3_general.pkl','rb'))
print("E3 keys:", list(E3.keys()))
print("nl count:", len(E3['nl']), " loc count:", len(E3['loc']))
print("first nl entries:")
for (key,co) in E3['nl'][:6]:
    print("  ", key, "->", co)
print("first loc entries:")
for (key,co) in E3['loc'][:6]:
    print("  ", key, "->", co)

Jtab = pickle.load(open('Jtab_t.pkl','rb'))
print("\nJtab entries:", len(Jtab))
for i,(k,v) in enumerate(Jtab.items()):
    print("  ", k, "->", v)
    if i>=5: break
