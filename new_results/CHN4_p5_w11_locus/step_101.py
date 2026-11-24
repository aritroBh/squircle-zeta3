import time, os
time.sleep(160)
print(open('gamma_final.log').read() if os.path.exists('gamma_final.log') else "not yet")
