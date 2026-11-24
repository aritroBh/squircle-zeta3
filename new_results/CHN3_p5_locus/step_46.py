import time
src = open('engine.py').read()
print(src[2000:])
time.sleep(20)
print("\n--- engine12.log ---")
print(open('engine12.log').read())
