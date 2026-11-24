import subprocess, os
base="../SRF3_octahedron_radius_proof"
out = subprocess.run(["grep","-rl","3*pi/2\|3pi/2\|phi_d\|mixing", base, "--include=step_*.py"], capture_output=True, text=True)
print(out.stdout)
out2 = subprocess.run(["grep","-rn","phi", base+"/progress.log"], capture_output=True, text=True)
print(out2.stdout[:2000])
