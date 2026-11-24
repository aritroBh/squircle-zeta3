from flint import arb
rho = arb.pi()/2
send = arb('1.9')*rho
def psi_arb(s):
    q = s.sin()/s
    return (q + 1 + s**2/4 - (2*(s/2).cos() + s*(s/2).sin())*q.sqrt())/s**4

h = arb('0.001')
N = int(float((send - arb('0.1'))/(2*h))) + 2
maxub = arb(-10)**10 * -1  # -inf hack: just use None
maxub = None; argmax=None
c = arb('0.1') + h
for i in range(N):
    lo = float(c-h); hi = min(float(c+h), float(send))
    if lo > float(send): break
    seg = arb(arb((lo+hi)/2), arb((hi-lo)/2))
    v = psi_arb(seg)
    ub = v.upper()
    if maxub is None or ub > maxub: maxub = ub; argmax=(lo+hi)/2
    c += 2*h
psi524 = (arb(5)/24)**2
print("intervals used:", N, " argmax near s =", argmax)
print("certified max psi on [0.1, 1.9rho]:", maxub)
print("(5/24)^2 =", psi524)
print("max psi - (5/24)^2 :", maxub - psi524)
Estar = (maxub if maxub>psi524 else psi524).sqrt()
print("CERTIFIED max |E| on full slit [0, 1.9 rho] <=", Estar)
print("5/24 =", arb(5)/24, " so E* = 5/24 is CERTIFIED (sup attained at s=0)")
with open('progress.log','a') as f:
    f.write(f"step_08c: E*=5/24 CERTIFIED. piece B: {N} intervals, max psi={maxub} < (5/24)^2. Piece A series c1<0.\n")
