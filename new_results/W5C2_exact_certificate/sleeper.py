import time
for i in range(12):
    time.sleep(10)
    open('sleeper.log','a').write(f'tick {i} {time.ctime()}\n')
