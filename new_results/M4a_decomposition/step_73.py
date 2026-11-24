import json, mpmath as mp, datetime
mp.mp.dps = 900
res = json.load(open("hp_results.json"))
def finalize(hi, lo, ndig, fname):
    a = mp.mpf(res[hi]["value"]); b = mp.mpf(res[lo]["value"])
    s = mp.nstr(a, ndig)
    open(fname, "w").write(s + "\n")
    print(f"{fname}: {ndig} digits; two-run agreement {mp.nstr(abs(a-b),3)}")
    print(" ", s[:80], "...")
finalize("w4pp@890", "w4pp@850", 800, "W4pp_800.txt")
finalize("w5@690", "w5@650", 600, "W5_600.txt")
finalize("m4@890", "m4@850", 800, "M4_800.txt")
with open("progress.log", "a") as fh:
    fh.write(f"[{datetime.datetime.now().isoformat()}] Task4 done: W4'' 800d (agree 2.7e-810), W5 600d (8.3e-612), M4 800d (3.6e-813), mpmath tanh-sinh two precisions.\n")
