import subprocess
out = subprocess.run(['grep','-rl','DIFF','<project>/research/08_theory/OCT6_octagon_closed'], capture_output=True, text=True).stdout
print(out)
out2 = subprocess.run(['grep','-l','A(x)','<project>/research/08_theory/OCT6_octagon_closed/*.py'], capture_output=True, text=True).stdout
print(out2)
