from flint import arb
def psi_arb(s):
    q = s.sin()/s
    return (q + 1 + s**2/4 - (2*(s/2).cos() + s*(s/2).sin())*q.sqrt())/s**4
for c in [0.101, 0.5, 1.0, 2.0, 2.98]:
    seg = arb(arb(c), arb(0.0005))
    print(c, psi_arb(seg))
