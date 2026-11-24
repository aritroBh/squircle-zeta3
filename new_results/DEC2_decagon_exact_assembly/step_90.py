import subprocess, time, os
p = "<project>/research/08_theory/PRV1_depth2_proofs"
out = subprocess.run(['head','-40', f'{p}/verify_full.py'], capture_output=True, text=True).stdout
print(out)
print("size:", os.path.getsize(f'{p}/verify_full.py'))
