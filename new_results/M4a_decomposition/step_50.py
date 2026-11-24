import os, hashlib
BASE = "<project>"
OUT = os.path.join(BASE, "paper/squircle_note/outputs")
for f in ["m3_tails_exact.pkl", "m4_tails_exact.pkl"]:
    h1 = hashlib.sha256(open(os.path.join(OUT, f), "rb").read()).hexdigest()[:16]
    h2 = hashlib.sha256(open(f, "rb").read()).hexdigest()[:16]
    print(f, h1, h2, "SAME" if h1 == h2 else "DIFFERENT")
