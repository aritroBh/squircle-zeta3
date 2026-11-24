base = "<project>/research/"
import pickle
wc = pickle.load(open(base + "08_theory/COR5_corner_classicality/word_classes.pkl", "rb"))
for k in sorted(wc.keys()):
    print(k, "naive", tuple(str(x) for x in wc[k][0]), "true", tuple(str(x) for x in wc[k][1]), "denom", wc[k][2])
