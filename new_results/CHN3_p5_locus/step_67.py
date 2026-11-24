import time, os
time.sleep(280)
print(open('eng12b.log').read())
print([f for f in ['eng12.pkl','vw12.pkl','forms1112.pkl'] if os.path.exists(f)])
try: os.kill(int(open('eng12b.pid').read()),0); print("alive")
except Exception as e: print(e)
