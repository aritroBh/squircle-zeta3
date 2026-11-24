import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = open("verify_oct5d3.py").read()
print("count of 'def main()':", src.count("def main()"))
print("count of \"__main__':\":", src.count("__main__':"))
print("count of 'ALL CHECKS PASS':", src.count("ALL CHECKS PASS"))
print("count of '== A. exact':", src.count("== A. exact"))
print("last 300 chars:", repr(src[-300:]))
