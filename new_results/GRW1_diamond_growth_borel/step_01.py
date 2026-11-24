import os
base = "<project>/"
for d in ["research/08_theory/GRW1_diamond_growth_borel",
          "research/08_theory/D11_rhombus_resurgence",
          "research/08_theory/D11b_stokes_numeric",
          "research/08_theory/D11b2_stokes_jump",
          "research/08_theory/D11b3_jump_fit",
          "research/08_theory/D11c_zeta_first_order",
          "research/08_theory/D13c_two_ends"]:
    p = os.path.join(base, d)
    print("="*80)
    print(d, "EXISTS" if os.path.isdir(p) else "MISSING")
    if os.path.isdir(p):
        for f in sorted(os.listdir(p)):
            print("   ", f)
