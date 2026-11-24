CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
for tag in ["mut_a", "mut_b", "mut_c", "mut_d"]:
    print(f"===== {tag} full =====")
    print(open(CWD + "/" + tag + ".out").read().split("--- STDERR ---")[0])
