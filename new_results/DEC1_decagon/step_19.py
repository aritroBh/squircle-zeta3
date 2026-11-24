import subprocess
d='research/08_theory/OCT6_octagon_closed/oct5/'
out=subprocess.run(['grep','-l','E3_octagon_exact',*[d+f'step_{i:02d}.py' for i in range(1,61)]],capture_output=True,text=True)
print(out.stdout, out.stderr)
out=subprocess.run(['grep','-rn','def build_E3\|E3_exact\|def E3',d],capture_output=True,text=True)
print(out.stdout[:3000])