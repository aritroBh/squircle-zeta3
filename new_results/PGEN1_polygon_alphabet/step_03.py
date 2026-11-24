base = "<project>"
import os
for d in ["research/08_theory/OCT9c_octagon_close", "research/08_theory/OCT9d_octagon_short",
          "research/08_theory/DOD5_dodecagon_termwise", "research/08_theory/DEC1_decagon",
          "research/08_theory/K1e_hexagon_closed", "research/08_theory/OCT8_symbol_test",
          "research/08_theory/OCT6_octagon_closed", "research/08_theory/OCT4_octagon_letters",
          "research/08_theory/CX1_octagon_depth2_proof", "research/08_theory/PRV1_depth2_proofs"]:
    p = os.path.join(base, d)
    print("=====", d)
    try:
        print(os.listdir(p))
    except Exception as e:
        print("ERR", e)
