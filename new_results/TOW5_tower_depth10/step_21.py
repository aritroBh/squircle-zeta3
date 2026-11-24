with open('<project>/research/08_theory/WTS_weight_k_structure/theory/run_kernel_k.py') as f:
    print(f.read())
print("="*80)
# How did TOW4 check k=5 exceptional? grep step files
import subprocess
out = subprocess.run(['grep','-l','kern_5','<project>/research/08_theory/TOW4_tower_depth8/step_7'+str(i)+'.py' for i in range(0,10)], capture_output=True, text=True)
print(out.stdout, out.stderr)
