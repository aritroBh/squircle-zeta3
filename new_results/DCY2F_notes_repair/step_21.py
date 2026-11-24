import re
base = "<project>/research/08_theory"
with open(f"{base}/DCY2R_referee_varpi1/transcript.json") as f:
    trr = f.read()
print("len:", len(trr))
for pat in ["res = abs(M - T1)", "res = abs(M + T1)", "res = abs(M - (-T1))", "abs(M - (-T1))", "target_Im"]:
    idxs = [m.start() for m in re.finditer(re.escape(pat), trr)]
    print(f"{pat!r}: {len(idxs)}")
# Show context where referee quotes DCY2's step_06 code
for m in list(re.finditer(re.escape("abs(M - (-T1))"), trr))[:3]:
    i=m.start()
    print("..."+trr[i-500:i+300]+"...")
    print("="*100)
