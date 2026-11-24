import json, pickle
import mpmath as mp
tails4 = pickle.load(open("m4_tails_exact.pkl", "rb"))   # BEFORE dps (import resets dps!)
mp.mp.dps = 950
W5 = mp.mpf(open("W5_600.txt").read().strip())
W4pp = mp.mpf(open("W4pp_800.txt").read().strip())
M4 = mp.mpf(open("M4_800.txt").read().strip())
W34 = tails4["W3p"].num()
T = tails4["T"].num()
s = W5 + W4pp + W34 + T
print("W5+W4''+W3^4+T - M4 (at ~600 shared digits) =", mp.nstr(s - M4, 6))

# Also verify the W4'' closed form against the independent G-expression value (real_g, ~200 digits)
import subprocess, sys
res = json.load(open("hp_results.json"))
chk = json.load(open("gpl_w4pp_check.json"))
print("W4'' G-expression (real_g) earlier check vs 90-digit quad:", chk["err_vs_quad"])
chk5 = json.load(open("gpl_w5_check.json"))
print("W5 G-expression (real_g) earlier check vs 90-digit quad:", chk5["err_vs_quad"])
wv4 = json.load(open("wv_w4pp_check.json")); wv5 = json.load(open("wv_w5_check.json"))
print("W4'' word-vector check:", wv4["err_vs_quad"], "| W5 word-vector check:", wv5["err_vs_quad"])
import datetime
with open("progress.log", "a") as fh:
    fh.write(f"[{datetime.datetime.now().isoformat()}] Task5 done: W4'' relation stable 700/1000 (sign), closed form agrees 1.07e-1060. Final sum check printed.\n")
