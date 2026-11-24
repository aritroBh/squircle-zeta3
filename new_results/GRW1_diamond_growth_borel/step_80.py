wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
src = open(wd+"/proof.md").read()
i = src.index("**Lemma 3.3 (square).**")
j = src.index("## Part III")
print(src[i:j])
