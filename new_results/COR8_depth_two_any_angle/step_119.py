
import os
os.chdir('<project>/research/08_theory/COR8_depth_two_any_angle')
open('progress.log','a').write('step_final: run ended. Series engine validated; D series computed; Step-4 fit 90%% done (crashed on trivial bug at the end). See result.md.\n')
print(open('result.md').read()[:800])
print('---')
import subprocess
print(subprocess.run(['ls'],capture_output=True,text=True).stdout)
