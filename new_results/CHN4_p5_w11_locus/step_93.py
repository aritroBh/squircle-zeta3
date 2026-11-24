import time, os
time.sleep(150)
print(open('gamma_w12_extra.log').read() if os.path.exists('gamma_w12_extra.log') else "not yet")
print("out:", open('run_gamma_w12_extra.out').read() if os.path.exists('run_gamma_w12_extra.out') else "none")
