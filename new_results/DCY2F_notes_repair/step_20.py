import re
base = "<project>/research/08_theory"
with open(f"{base}/DCY2_decay_first_correction/transcript.json") as f:
    tr = f.read()
for m in list(re.finditer(re.escape("abs(M - T1)"), tr))[:3]:
    i = m.start()
    print("..."+tr[i-700:i+150]+"...")
    print("="*100)
