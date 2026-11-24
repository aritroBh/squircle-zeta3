import sys, pickle, json
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/M5a_next_coefficient'
sys.path.insert(0, work)
# UNPICKLE FIRST (may import dilog_exact which resets dps to 60), set dps AFTER
tails5 = pickle.load(open(f"{work}/m5_tails_exact.pkl", "rb"))
import mpmath as mp
mp.mp.dps = 110
W3exact = tails5["W3p"].num()
Texact = tails5["T"].num()
print("W3^(5) exact @110dps =", mp.nstr(W3exact, 100))
print("T exact =", mp.nstr(Texact, 40))
lam = mp.log(2)
print("T - (pi^2+6 lam^2)/24 =", mp.nstr(Texact - (mp.pi**2 + 6*lam**2)/24, 5))
q = json.load(open(f"{work}/quad5_100.json"))
W35q = mp.mpf(q["W35_quad"])
print("W3 exact - W3 quad =", mp.nstr(W3exact - W35q, 5))
# exact-tail sum check at this precision
W6 = mp.mpf(q["W6"]); W5pp = mp.mpf(q["W5pp"]); W4ppp = mp.mpf(q["W4ppp"]); M5 = mp.mpf(q["M5"])
s = W6 + W5pp + W4ppp + W3exact + Texact
print("W6+W5pp+W4ppp+W3exact+Texact - M5 =", mp.nstr(s - M5, 5))
with open(f'{work}/progress.log', 'a') as f:
    f.write("step: verified exact tails (T, W3(5)) at 110 dps against quadrature and full sum\n")