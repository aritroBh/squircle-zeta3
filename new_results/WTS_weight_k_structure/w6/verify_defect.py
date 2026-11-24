"""verify step: k=5 defect delta = W_6(5) - p(5) = -Delta^5 over k=5..10, from OWN 110-digit values; compare with the quoted
-0.219366913442151079860842560263653086600349859821953115922874 and with the claim 'Delta^5 on k=5..10 = 0.2193669134421510798608426'.
usage: verify_defect.py OWNDIR"""
import sys, os, json
sys.dont_write_bytecode = True
import mpmath as mp
OWN = sys.argv[1]
mp.mp.dps = 130
print("mp.mp.dps =", mp.mp.dps)
y = [mp.mpf(json.load(open(os.path.join(OWN, f"own_k{k}_d5_110.json")))["value"]) for k in range(5, 11)]
d5 = sum((-1) ** (5 - m) * mp.binomial(5, m) * y[m] for m in range(6))
print("Delta^5 (k=5..10) own  =", mp.nstr(d5, 60))
print("delta = W_6(5)-p(5) own =", mp.nstr(-d5, 60))
q = mp.mpf("-0.219366913442151079860842560263653086600349859821953115922874")
print("quoted delta - own delta =", mp.nstr(q + d5, 4))
