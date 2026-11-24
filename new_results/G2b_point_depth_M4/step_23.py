
base = "<project>/research/08_theory"
import subprocess
for f in ["hex_W4.txt", "hex_F4.txt", "hex_W3.txt", "hex_W2.txt", "hex_Psi1.txt", "hex_Psi2.txt", "hex_R3_2.txt", "hex_R3_3.txt"]:
    s = open(base + "/COR5_corner_classicality/" + f).read().strip()
    print(f, "->", s[:60], "...", "len:", len(s))
