import time, os, pickle
while not os.path.exists("gsym_tab.pkl"):
    time.sleep(15)
print(open("gsym_build.log").read()[-200:])
tab = pickle.load(open("gsym_tab.pkl","rb"))
print("words:", len(tab))
