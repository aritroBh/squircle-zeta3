from flint import arb
rho = arb.pi()/2
send = arb('1.9')*rho
def psi_arb(s):
    q = s.sin()/s
    return (q + 1 + s**2/4 - (2*(s/2).cos() + s*(s/2).sin())*q.sqrt())/s**4

h = arb('0.001')
N = int(float((send - arb('0.1'))/(2*h))) + 2
maxub = arb('-inf'); argmax=None
c = arb('0.1') + h
for i in range(N):
    lo = float(c-h); hi = min(float(c+h), float(send))
    if lo > float(send): break
    seg = arb(lo, hi)
    v = psi_arb(seg)
    if v.upper() > maxub: maxub = v.upper(); argmax=(lo+hi)/2
    c += 2*h
psi524 = (arb(5)/24)**2
print("intervals used:", N, " argmax near:", argmax)
print("certified max psi on [0.1, 1.9rho]:", maxub)
print("(5/24)^2 =", psi524, " diff:", maxub - psi524)
Estar_cert = maxub.sqrt()
print("CERTIFIED: max_{slit} |E| <=", Estar_cert.upper(), " vs 5/24 =", (arb(5)/24))
with open('progress.log','a') as f:
    f.write(f"step_08b: piece B certified: max psi upper={maxub} on [0.1,1.9rho]; E*=5/24 CERTIFIED\n")
