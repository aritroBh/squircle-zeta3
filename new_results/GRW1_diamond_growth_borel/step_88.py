wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
src = open(wd+"/proof.md").read()
i = src.index("**Theorem 2 (Borel summability — PROVED).**")
j = src.index("## Part V")
print(src[i:j])
