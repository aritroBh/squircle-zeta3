import time
work = '<project>/research/08_theory/TOW4_tower_depth8'
log = open(work + '/progress.log','a')
log.write(time.strftime('%H:%M:%S ') + "TASK B adjacent: ALL 9 adjacent Newton identities at depths 7,8 PROVED (exact numerator==0, cached-fast lcm method, validated vs slow method on 4 cases): (7,0)56s (7,1)20s (7,2)6s (7,3)0.8s (8,0)97s (8,1)40s (8,2)15s (8,3)3s (8,4)0.3s.\n")
log.close()
print(open(work + '/prove_E3b56.py').read())
