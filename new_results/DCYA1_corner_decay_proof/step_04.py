import subprocess
# Read the obs:Mk section (1243-1290)
print("=== Section: How Fast Do the M_k Decay? (1243-1290) ===")
print(subprocess.run(['sed', '-n', '1243,1295p', '<project>/paper/ieee/main.tex'], capture_output=True, text=True).stdout)
