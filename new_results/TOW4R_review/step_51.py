
import os, re
base = "<project>"
txt = open(base + "/research/08_theory/WTS_weight_k_structure/theory/result.md").read()
for m in re.finditer(r"\*\*Lemma 2", txt):
    print(txt[m.start():m.start()+1500])
    print("="*80)
