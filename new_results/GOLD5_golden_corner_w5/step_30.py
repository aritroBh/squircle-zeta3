import os, subprocess
print(sorted(os.listdir('../M4e_w5_lchi')))
out = subprocess.run(['head', '-c', '4000', '../M4e_w5_lchi/result.md'], capture_output=True, text=True)
print(out.stdout)
