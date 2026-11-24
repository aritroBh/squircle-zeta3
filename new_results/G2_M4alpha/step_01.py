import os

# Explore the working directory and the input directories
base = "<project>/research/"
for d in ["10_frontier/G1_zeta3_any_angle", "10_frontier/W1_corner_family/B_taylor_pi4", "08_theory/M4c_w5_closed_form", "08_theory/M4e_w5_lchi", "08_theory/COR6_true_point_symbols"]:
    p = os.path.join(base, d)
    print("="*80)
    print(p)
    try:
        for f in sorted(os.listdir(p)):
            print("   ", f)
    except Exception as e:
        print("ERR", e)
