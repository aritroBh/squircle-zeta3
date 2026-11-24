import time
time.sleep(1200)
s = open("val96_out.txt").read()
print("sequential done?", "sequential chamber_S done" in s)
print(s[-600:])
