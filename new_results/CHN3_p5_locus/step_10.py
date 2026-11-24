print("--- planner_p5_fft.py ---")
print(open('../CHN2_chain_planes/planner/planner_p5_fft.py').read()[:3000])
print("\n--- my existing step files (quick peek) ---")
import os
for f in sorted(f for f in os.listdir('.') if f.startswith('step_')):
    print(f"### {f} ({os.path.getsize(f)}b)")
    print(open(f).read()[:400])
    print()
