base='<project>/research/08_theory'
print(open(base+'/COR5_corner_classicality/F4_relation.json').read())
import subprocess
print(subprocess.run(['grep','-n','Psi2\|W4 =\|W4=\|zeta(3)\|M3 =\|M3=','-i',base+'/COR5_corner_classicality/result.md'],capture_output=True,text=True).stdout[:2500])
