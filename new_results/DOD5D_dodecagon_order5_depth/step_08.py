base = "<project>/research/08_theory/"
for f in ["M4c_w5_closed_form/w5_symbol.py", "M4c_w5_closed_form/w5_cobracket.py", "M4c_w5_closed_form/planner/planner_symbol_pairing.py"]:
    print("="*100); print(f); print("="*100)
    with open(base+f) as fh:
        print(fh.read())
