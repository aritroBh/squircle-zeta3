
import os, re
base = "<project>"
txt = open(base + "/research/08_theory/WTS_weight_k_structure/theory/result.md").read()
i = txt.find("Lemma 3")
# print full sections around lemma statements in the body (not table)
for m in re.finditer(r"\*\*Lemma [23]", txt):
    print("="*80)
    print(txt[m.start():m.start()+2200])
