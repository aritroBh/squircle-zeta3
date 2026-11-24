d = "research/08_theory/K4b_cubo_finish"
with open(d + "/k4_quad.py") as f:
    print(f.read())
print("="*80)
for s in ["step_01.py","step_02.py","step_03.py"]:
    print("###", s)
    with open(d + "/" + s) as f:
        print(f.read()[:6000])
    print("-"*80)
