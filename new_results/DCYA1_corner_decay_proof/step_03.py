import subprocess
# Find the corner section in the paper
out = subprocess.run(['grep', '-n', 'kappa_1', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout
print("=== kappa_1 mentions ===")
print(out)
out2 = subprocess.run(['grep', '-n', 'Corners of Any Angle\|thm:Malpha\|app:corner\|obs:Mk\|How Fast Do the\|M_k(\\\\alpha)\|M_k(', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout
print("=== section/theorem mentions ===")
print(out2)
