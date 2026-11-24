src = open("polylog_exact.py").read()
import re
i = src.find("def solve_li3_linear")
print(src[i:i+2600])
print("==== context of line 303 (workdps 60) ====")
lines = src.split("\n")
print("\n".join(lines[270:330]))
