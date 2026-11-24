import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
d1p = pickle.load(open('dod5_task1.pkl','rb'))
TERMS, vals = d1p['TERMS'], d1p['vals']
mp.mp.dps = 620
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
rv = pickle.load(open(D3+"R_all12_vals620.pkl","rb"))
TRIPLES = [(1,2,5),(1,2,6),(1,5,6),(2,1,5),(2,1,6),(2,5,6),(5,1,2),(5,1,6),(5,2,6),(6,1,2),(6,1,5),(6,2,5)]
print("triple:  |Rv(sym) - Rq_accurate|   (was: |Rv - Rq_gdswaps|)")
for (i,j,k) in TRIPLES:
    Rq_acc = atoms['T3a'][(i,j,k)] - (vals[('SwG',(i,j,k))]-vals[('SwD',(i,j,k))])/atoms['X'][i]
    dnew = abs(rv['Rv'][(i,j,k)] - Rq_acc)
    dold = abs(rv['Rv'][(i,j,k)] - rv['Rq'][(i,j,k)])
    print("  %-10s %s   (was %s)"%(str((i,j,k)), mp.nstr(dnew,3), mp.nstr(dold,3)))
