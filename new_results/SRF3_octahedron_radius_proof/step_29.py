
import subprocess, os, time
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
# check engine progress
print(open(wd+"/engine_pid.txt").read())
print(subprocess.run(["ps","-p",open(wd+"/engine_pid.txt").read().strip(),"-o","etime,stat"],capture_output=True,text=True).stdout)
print(open(wd+"/engine_K120.log").read()[-300:])
