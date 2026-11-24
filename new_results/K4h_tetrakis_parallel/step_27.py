s = open("val96_out.txt").read()
import re
for line in s.splitlines():
    if "sequential" in line or "parallel done" in line:
        print(line)
