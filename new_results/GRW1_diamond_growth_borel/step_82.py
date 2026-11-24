wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
src = open(wd+"/proof.md").read()
i = src.index("*Singular directions:*")
j = src.index("---", i)
print(src[i:j])
