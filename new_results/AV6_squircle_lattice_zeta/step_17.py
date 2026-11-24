src = open("squircle.py").read()
print(".real" in src, src.count(".real"))
# inspect the Jint function
import re
m = re.search(r"def Jint.*?(?=\ndef )", src, re.S)
print(m.group(0))
