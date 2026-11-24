import subprocess
base="../SRF3_octahedron_radius_proof"
out = subprocess.run(["grep","-n","phi\|3pi/2\|3\*pi/2\|mixing\|e\*\*(-z\|exp(-z", base+"/step_48.py"], capture_output=True, text=True)
print(out.stdout[:5000])
