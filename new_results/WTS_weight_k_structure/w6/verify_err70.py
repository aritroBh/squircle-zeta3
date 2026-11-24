"""verify step: true error of the 70-dps run (vals/) and of the 110-dps run (vals_110/) measured against OWN 110-digit values
(own quadrature, mp.mp.dps=120).  usage: verify_err70.py OWNDIR"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
OWN = sys.argv[1]
mp.mp.dps = 130
print("mp.mp.dps =", mp.mp.dps)
for k in range(5, 16):
    fo = os.path.join(OWN, f"own_k{k}_d5_110.json")
    if not os.path.exists(fo):
        continue
    o = mp.mpf(json.load(open(fo))["value"])
    s70 = mp.mpf(json.load(open(os.path.join(W, "vals", f"k{k}_d5.json")))["value"])
    s110 = mp.mpf(json.load(open(os.path.join(W, "vals_110", f"k{k}_d5.json")))["value"])
    e70 = json.load(open(os.path.join(W, "vals", f"k{k}_d5.json")))["err_est"]
    print(f"k={k}: |theirs70 - own| = {mp.nstr(abs(s70 - o), 3)} (their err_est {e70});  |theirs110 - own| = {mp.nstr(abs(s110 - o), 3)}")
