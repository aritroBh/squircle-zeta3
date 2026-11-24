import time, os, pickle
while not os.path.exists("gsym_tab.pkl"):
    time.sleep(10)
print(open("gsym_build.log").read()[-300:])
