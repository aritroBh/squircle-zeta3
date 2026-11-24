g1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
for tag in ["pi_2", "2pi_3", "pi_3"]:
    with open(g1 + f"values_{tag}.txt") as f:
        print(f"--- values_{tag}.txt ---")
        print(f.read())
