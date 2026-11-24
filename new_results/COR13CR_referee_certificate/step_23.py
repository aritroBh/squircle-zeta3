import os
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
print(open(G1 + "/flat_symbolic_proof.out").read()[:3000] if os.path.exists(G1+"/flat_symbolic_proof.out") else "no file")
