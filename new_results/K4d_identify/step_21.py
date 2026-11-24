import os, time
WORK = 'research/08_theory/K4d_identify'
os.chdir(WORK)
print("out exists:", os.path.exists('_sleep_out.txt'))
if os.path.exists('_sleep_out.txt'):
    print("content:", open('_sleep_out.txt').read())
print("log:", open('_sleep_log.txt').read() if os.path.exists('_sleep_log.txt') else 'NONE')
