import subprocess
out=subprocess.run(['grep','-rn','-e','38','-e','84','../M4e_w5_lchi/result_part2.md'],capture_output=True,text=True).stdout
print(out[:3000])
