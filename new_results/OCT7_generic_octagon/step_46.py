src = open("mixed_engine_gen.py").read()
old_d5 = """D5 = sp.expand((sp.log((k1-k2)/k1)*ell1 + L2mh - L2m2t)/k1)   # log(2/3) -> expanded below
D5 = sp.expand(D5.subs(sp.log(sp.Rational(2,3)), l3-l2))"""
new_d5 = """log_23 = l3 - l2     # log(2/3) = log((k1-k2)/k1), as symbols
D5 = sp.expand((log_23*ell1 + L2mh - L2m2t)/k1)"""
assert old_d5 in src; src = src.replace(old_d5, new_d5)

old_d7 = """D7 = (ell1**2*sp.log(1-u0/R1) - 2*ell1*L2m3h - 2*L3m3h + ell1**3/3 + 2*L3m2)/k2
D7 = sp.expand(D7.subs(sp.log(1-u0/R1), l5-l2))               # log(5/2)"""
new_d7 = """log_52 = l5 - l2     # log(5/2) = log(1-u0/R1)
D7 = sp.expand((ell1**2*log_52 - 2*ell1*L2m3h - 2*L3m3h + ell1**3/3 + 2*L3m2)/k2)"""
assert old_d7 in src; src = src.replace(old_d7, new_d7)

old_d8 = """D8 = (ell2**2*sp.log(1-u0b/R2) - 2*ell2*L2_3f - 2*L3_3f + ell2**3/3 + 2*L3_2t)/k1
D8 = sp.expand(D8.subs(sp.log(1-u0b/R2), l2-l5))              # log(2/5)"""
new_d8 = """log_25 = l2 - l5     # log(2/5) = log(1-u0b/R2)
D8 = sp.expand((ell2**2*log_25 - 2*ell2*L2_3f - 2*L3_3f + ell2**3/3 + 2*L3_2t)/k1)"""
assert old_d8 in src; src = src.replace(old_d8, new_d8)
open("mixed_engine_gen.py","w").write(src)
print("patched. verifying no bare log atoms remain in D constants:")
import importlib, mixed_engine_gen as E
importlib.reload(E)
import sympy as sp
ok = True
for name in ['D1','D2','D3','D4','D5','D6','D7','D8','D9','D10']:
    D = getattr(E, name)
    bad = [a for a in sp.preorder_traversal(D) if isinstance(a, sp.log)]
    if bad: print(name, "HAS log atoms:", bad); ok=False
print("clean" if ok else "NOT clean")
