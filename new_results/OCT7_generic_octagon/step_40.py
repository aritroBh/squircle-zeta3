src = open("mixed_engine_gen.py").read()
old = "    return sp.expand((-ell1**2*ell2 - 2*k1*E1 - k2*E2)/(n+1))"
new = "    return sp.expand((+ell1**2*ell2 - 2*k1*E1 - k2*E2)/(n+1))"
assert old in src
src = src.replace(old, new)
open("mixed_engine_gen.py","w").write(src)
print("fixed M3 boundary sign")

# also fix the Jloc test case in step_04 (0,2,2) -> in-range
s4 = open("step_04.py").read()
s4 = s4.replace("(3,3,0),(0,2,2),(3,1,2)", "(3,3,0),(1,2,1),(3,1,2)")
open("step_04.py","w").write(s4)
print("fixed step_04 test set")
