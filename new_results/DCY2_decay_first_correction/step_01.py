import os
base = "<project>/"
# List working dir and relevant dirs
for d in ["research/08_theory/DCY2_decay_first_correction",
          "research/08_theory/CX6_decay_rest",
          "research/08_theory/D13s_fix_recertify",
          "research/08_theory/D13u_close_h3",
          "research/08_theory/D13f_sweep"]:
    p = os.path.join(base, d)
    print("="*20, d, "exists:", os.path.exists(p))
    if os.path.exists(p):
        for f in sorted(os.listdir(p)):
            print("   ", f, os.path.getsize(os.path.join(p,f)))
