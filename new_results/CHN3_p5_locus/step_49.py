import time
time.sleep(240)
print(open('engine12.log').read())
import os
print("forms12.pkl:", os.path.exists('forms12.pkl'), os.path.getsize('forms12.pkl') if os.path.exists('forms12.pkl') else '')
