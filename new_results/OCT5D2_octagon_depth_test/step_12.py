import os
p = "<project>/paper/squircle_note/scripts/"
print(os.path.exists(p), os.listdir(p) if os.path.exists(p) else "")
print(open(p+"qsqrt2.py").read()[:6000])
