import subprocess, sys, json
r = subprocess.run([sys.executable, "step_hp.py", "m4", "850"], capture_output=True, text=True, timeout=2700)
print(r.stdout); print("ERR", r.stderr[-300:])

# Agreement checks
import mpmath as mp
mp.mp.dps = 900
res = json.load(open("hp_results.json"))
def agree(hi, lo, ndig):
    a = mp.mpf(res[hi]["value"]); b = mp.mpf(res[lo]["value"])
    d = abs(a-b)
    print(f"{hi} vs {lo}: |diff| = {mp.nstr(d, 3)}  (target {ndig} digits)")
    return d
agree("w4pp@890", "w4pp@850", 800)
agree("w5@690", "w5@650", 600)
agree("m4@890", "m4@850", 800)
