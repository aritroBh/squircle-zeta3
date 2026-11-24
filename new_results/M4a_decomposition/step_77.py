import subprocess, sys, json
r = subprocess.run([sys.executable, "step_pslq2.py", "W4pp_1060.txt", "1000"], capture_output=True, text=True, timeout=2700)
print(r.stdout[-3300:])
print("ERR", r.stderr[-500:])
# stability comparison
j7 = json.load(open("pslq2_w4pp_700.json"))
j10 = json.load(open("pslq2_w4pp_1000.json"))
same = j7["relation"] == j10["relation"]
print("RELATION IDENTICAL at 700 and 1000 digits:", same)
if not same:
    for k in set(j7["relation"]) | set(j10["relation"]):
        a, b = j7["relation"].get(k, 0), j10["relation"].get(k, 0)
        if a != b:
            print("  DIFF", k, a, b)
