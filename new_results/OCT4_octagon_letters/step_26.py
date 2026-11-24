import os, sys, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 300
r2 = mp.sqrt(2); I = 1j
# 8 real cross-ratio reps
real_pts = {'Li3(-1)': -1, 'Li3(-r2/2)': -r2/2, 'Li3(1-r2)': 1-r2, 'Li3(1-r2/2)': 1-r2/2,
            'Li3(r2-1)': r2-1, 'Li3(1/2)': mp.mpf('0.5'), 'Li3(2-r2)': 2-r2, 'Li3(r2/2)': r2/2}
# 12 canonical complex letters (Im>0,|z|<=1)
def C(re,im): return mp.mpc(re,im)
s8 = mp.sin(mp.pi/8); c8=mp.cos(mp.pi/8)
complex_pts = {
 'CL1': C(*( (1-1j)*(-r2+1+1j)/4 ).real and (( (1-1j)*(-r2+1+1j)/4 ).real, ((1-1j)*(-r2+1+1j)/4).imag)),
}
# build them directly from sympy-derived forms
import sympy as sp
sr2=sp.sqrt(2); II=sp.I
forms = {
 'CL1': (1-II)*(-sr2+1+II)/4,
 'CL2': (1+II+sr2*II)/(3+II+sr2*(2+II)),
 'CL3': (-1+II+sr2*II)/(3+II+sr2*(2+II)),
 'CL4': sp.Rational(1,2)-II/2+sr2*II/2,
 'CL5': -sr2/2+sp.Rational(1,2)+II/2,
 'CL6': sp.Rational(1,2)+II/2,
 'CL7': 2*(1+sr2)/(2*sr2+3-II*(1+sr2)),
 'CL8': 2/(1-II*(1+sr2)),
 'CL9': sr2/4+sp.Rational(1,2)+sr2*II/4,
 'CL10': (1+sr2+II)/(1+sr2-II),
 'CL11': II,
 'CL12': (-sr2-1+II)/(1+sr2+II),
}
cvals={}
for k,v in forms.items():
    nv = complex(sp.N(v,40))
    cvals[k]=mp.mpc(nv.real,nv.imag)
    print(k, "-> |z|=",round(abs(cvals[k]),6)," arg/pi=",round(float(mp.arg(cvals[k])/mp.pi),4))
# save for reuse
import pickle
pickle.dump({'real_pts':{k:mp.nstr(v,50) for k,v in real_pts.items()},
             'complex_pts':{k:(mp.nstr(v.real,50),mp.nstr(v.imag,50)) for k,v in cvals.items()}},
            open('letter_pts.pkl','wb'))
log("step: defined 8 real + 12 complex letter points")
