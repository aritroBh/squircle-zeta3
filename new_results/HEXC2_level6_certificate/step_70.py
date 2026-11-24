base='<project>/research/08_theory'
import subprocess
out = subprocess.run(['grep','-n','Psi2\|W4\|zeta(3)\|M3 =','-i',base+'/COR5_corner_classicality/result.md'],capture_output=True,text=True).stdout
print(out[:3000])
