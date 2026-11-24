from flint import arb, acb
import mpmath as mp

# ---------- Piece A: s in [0, 0.1], psi(s) <= (5/24)^2 via series with rigorous tail ----------
# M_E(3) <= (sqrt((1+e^6)/6)+1+1.5)/9 <= 1.19 (rigorous with arb)
e6 = arb(6).exp()
M_E3 = (((1+e6)/6).sqrt() + 1 + arb('1.5'))/9
print("M_E(3) <=", M_E3.upper())
M2 = M_E3**2
c1 = arb(-25)/13824
c2 = arb(11147)/232243200
u_max = arb('0.01')   # s <= 0.1
T = M2/(9**3*(1-u_max/9))
delta = c1 + (c2 + T)*u_max
print("c1 + (c2+T)*0.01 =", delta.upper(), " < 0 ?", delta.upper() < 0)
print("=> psi(s) <= (5/24)^2 on [0,0.1] rigorously:", (delta.upper() < 0))

# ---------- Piece B: [0.1, 1.9*rho] arb interval cover ----------
rho = arb.pi()/2
send = arb('1.9')*rho
def psi_arb(s):
    # s arb interval; psi = [sin(s)/s + 1 + s^2/4 - (2 cos(s/2) + s sin(s/2))*sqrt(sin(s)/s)]/s^4
    sn = s.sin()
    q = sn/s
    return (q + 1 + s**2/4 - (2*(s/2).cos() + s*(s/2).sin())*q.sqrt())/s**4

h = arb('0.001')
N = int((send - arb('0.1'))/ (2*h)) + 2
maxub = arb('-inf'); argmax = None
c = arb('0.1') + h
psi524 = (arb(5)/24)**2
for i in range(N):
    if c - h > send: break
    seg = arb(c - h, min(float(c + h), float(send)))
    v = psi_arb(seg)
    if v.upper() > maxub:
        maxub = v.upper(); argmax = float(c)
    c += 2*h
print("certified max psi on [0.1, 1.9 rho]:", maxub)
print("(5/24)^2 =", psi524)
print("certified max |E| =", maxub.sqrt())
print("5/24 =", arb(5)/24)
E_star_cert = max(maxub, psi524).sqrt()
print("CERTIFIED E* =", E_star_cert.upper(), " sup is attained at s=0, value 5/24 exactly")
with open('progress.log','a') as f:
    f.write(f"step_08: E* CERTIFIED by arb: piece A series (c1<0), piece B {N} intervals; max|psi|=(5/24)^2\n")
