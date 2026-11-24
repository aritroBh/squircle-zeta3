wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
src = open(wd+"/proof.md").read()
i = src.index("## Lemma ledger")
print(src[i:])
import os
print("\n--- deliverables present:", [f for f in ["proof.md","result.md","result_nofinal.md","verify_grw1.py","progress.log"] if os.path.exists(os.path.join(wd,f))])
