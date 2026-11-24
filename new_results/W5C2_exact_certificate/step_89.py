import subprocess
print(subprocess.run(['bash','-c','grep -E "step_01b mp: feed done|step_01b mp: \|S\||step_01b mp: plateau" step_01b.out'], capture_output=True, text=True).stdout)
print(subprocess.run(['bash','-c','tail -3 step_01b.out'], capture_output=True, text=True).stdout)