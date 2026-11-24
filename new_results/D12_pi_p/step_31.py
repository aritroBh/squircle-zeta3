
import subprocess, time
d = "<project>/research/08_theory/D12_pi_p"
for f in ['step_04.py','step_02.py','step_03.py']:
    t=time.time()
    r = subprocess.run(['python3', f'{d}/{f}'], capture_output=True, text=True, timeout=900)
    print("="*20, f, f"({time.time()-t:.0f}s)")
    print(r.stdout[-1200:])
    if r.returncode!=0: print("STDERR:", r.stderr[-800:])
