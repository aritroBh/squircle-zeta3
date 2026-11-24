
for f in ["hex_Psi2.txt", "hex_Psi1.txt", "hex_W4.txt", "hex_W3.txt", "hex_W2.txt", "hex_F4.txt", "hex_R3_3.txt", "hex_R3_2.txt"]:
    txt = open(f"../COR5_corner_classicality/{f}").read().strip()
    print(f, ":", txt[:80], "... (%d chars)" % len(txt))
print()
m3 = open("../COR5_corner_classicality/M3_2pi3_1100.txt").read().strip()
print("M3_2pi3_1100.txt:", m3[:100], "... (%d chars)" % len(m3))
