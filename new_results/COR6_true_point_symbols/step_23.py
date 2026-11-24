base = "<project>/research/"
import pickle
wc = pickle.load(open(base + "08_theory/COR5_corner_classicality/word_classes.pkl", "rb"))
print(type(wc))
if isinstance(wc, dict):
    ks = list(wc.keys())[:6]
    for k in ks: print(k, "->", wc[k])
    print("n =", len(wc))
with open("step_01.py") as f:
    print(f.read()[:2000])
