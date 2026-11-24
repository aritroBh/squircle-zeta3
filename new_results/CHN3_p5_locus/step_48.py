import os, time
try:
    os.kill(91501, 0)
    print("engine12 process alive")
except Exception as e:
    print("process check:", e)
time.sleep(150)
print(open('engine12.log').read())
print("forms12.pkl:", os.path.exists('forms12.pkl'))
