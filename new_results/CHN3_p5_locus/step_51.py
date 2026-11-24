import time
time.sleep(280)
print(open('engine12.log').read())
import os
print("forms12.pkl:", os.path.exists('forms12.pkl'))
try:
    os.kill(91501,0); print("alive")
except Exception as e: print("dead?", e)
