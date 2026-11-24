from flint import arb
t = arb(0.5, 0.501)
print("t =", t)
print("t.sin() =", t.sin())
print("arb.sin(t) =", arb.sin(t))
s2 = arb(2)/3
print("s2=", s2, "s2.sin() =", s2.sin())
# maybe issue: arb(0.5,0.501) with floats? try mid/rad
t2 = arb(arb(0.5), arb(0.001))
print("t2 =", t2, "sin:", t2.sin())
