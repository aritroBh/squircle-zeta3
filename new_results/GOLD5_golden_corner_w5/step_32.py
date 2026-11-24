import subprocess
out = subprocess.run(['grep', '-n', 'class\|P5\|Li5 span\|li5_span\|classical class', '../G2b_point_depth_M4/step_17.py'], capture_output=True, text=True)
print(out.stdout[:1500])
print("=========")
src = open('../G2b_point_depth_M4/step_17.py').read()
print(src[:5000])
