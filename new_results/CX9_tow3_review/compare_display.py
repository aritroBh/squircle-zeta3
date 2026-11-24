"""Check the two displayed rho numerators, independently transcribed from proof_depth6.md."""
from pathlib import Path
import sys
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'TOW3_tower_depth6'))
import tow3lib as T
d,m,i=T.d,T.m,s.I
N5=(-16*i*d**7+304*i*d**6+(-80-2520*i)*d**5+160*d**4*m+(1040+11840*i)*d**4-2080*d**3*m+(-5140-34069*i)*d**3+(10280-240*i)*d**2*m+(11440+60081*i)*d**2+240*i*d*m**2+(-22880+960*i)*d*m+(-9600-59920*i)*d-960*i*m**2+19200*m+25200*i)
N6=(-64*d**9+1728*d**8+(-20976+480*i)*d**7-960*i*d**6*m+(150432-10080*i)*d**6+20160*i*d**5*m+(-702876+89040*i)*d**5+(-2880-178080*i)*d**4*m+(2221452-423840*i)*d**4+2880*d**3*m**2+(34560+847680*i)*d**3*m+(-4760069+1146870*i)*d**3-34560*d**2*m**2+(-138960-2293260*i)*d**2*m+(6677868-1670730*i)*d**2+(138960-1440*i)*d*m**2+(186480+3341460*i)*d*m+(-5524995+1019460*i)*d+960*i*m**3-186480*m**2-2038920*i*m+1984500)
assert s.cancel(T.arho(5)-N5/(480*T.Pi(2)))==0
assert s.cancel(T.arho(6)-N6/(5760*T.Pi(3)))==0
print('PASS displayed a5,a6 numerators equal submitted arho symbolically')
