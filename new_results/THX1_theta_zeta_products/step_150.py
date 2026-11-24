
import pickle
d89 = pickle.load(open('I_reductions_w89.pkl','rb'))
reds8 = d89['reds']; basis8 = d89['basis8']
for (m, j), coeffs in sorted(reds8.items()):
    terms = []
    for c, nm in zip(coeffs, basis8):
        if c:
            terms.append(f"{c}*{nm}" if c != 1 else nm)
    print(f"I_{{{m},{j}}} = " + " + ".join(terms))
