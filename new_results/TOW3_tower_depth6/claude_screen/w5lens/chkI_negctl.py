# negative control: perturb the proved e(6,1) at key (3,2) (and separately a_6) and confirm the chain checks detect it
import re,sys,os
sys.dont_write_bytecode=True
src=open('chkC_chain.py').read()
# perturb e(6,1) of key (3,2) by 1/1000 in the real part
src=src.replace("e00={key:QSI(*H[key][0][5]) for key in KEYS}",
 "E[((3,2),6,1)]=(E[((3,2),6,1)][0]+F(1,1000),E[((3,2),6,1)][1])\ne00={key:QSI(*H[key][0][5]) for key in KEYS}")
src=src.replace("for k in range(5,61):\n    KC[k]","for k in range(5,12):\n    KC[k]").replace("for k in range(5,22):\n    fn=","for k in range(5,12):\n    fn=")
src=src.replace("for k in range(5,61):\n        for key","for k in range(5,12):\n        for key") if False else src
exec(compile(src,'perturbed','exec'))
