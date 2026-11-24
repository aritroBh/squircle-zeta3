base = "<project>/research/08_theory/Q1S_vertex_symbol/"
for f in ["result.md", "q1s_stage1.py", "q1s_exact.py", "q1s_symbol.py"]:
    print("="*80)
    print("FILE:", f)
    print("="*80)
    print(open(base+f).read())
