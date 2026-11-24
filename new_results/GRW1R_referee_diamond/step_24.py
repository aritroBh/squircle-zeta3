from flint import arb
def psi_arb(s):
    q = s.sin()/s
    return (q + 1 + s**2/4 - (2*(s/2).cos() + s*(s/2).sin())*q.sqrt())/s**4
t = arb(0.5, 0.501)
print("q:", t.sin()/t)
v = psi_arb(t)
print("psi:", v)
print(type(v), v.upper() if hasattr(v,'upper') else 'no upper')
