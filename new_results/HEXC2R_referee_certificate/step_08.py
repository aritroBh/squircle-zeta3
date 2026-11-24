
src = open("../HEXC2_level6_certificate/hexc2_common.py").read().splitlines()
print(len(src))
print("\n".join(f"{i+1:4d} {l}" for i, l in enumerate(src[95:230])))
