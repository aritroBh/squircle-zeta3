import time
time.sleep(900)
s = open("val96_out.txt").read()
print(s[-800:])
print("---- sequential done?" , "sequential chamber_S done" in s)
