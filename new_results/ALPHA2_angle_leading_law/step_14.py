import pickle
from fractions import Fraction
with open('beta_tables.pkl','rb') as f: bt = pickle.load(f)
cvals = {'c12':Fraction(1,2),'c35':Fraction(3,5),'c13':Fraction(1,3),'c45':Fraction(4,5)}
for cname,c in cvals.items():
    r = (1-c)/(1+c)
    tab = bt[cname]
    print(f"===== {cname}: c={c}, r={r} =====")
    ds = sorted(set(k[0] for k in tab))
    for d in ds:
        print(f"--- d={d} ---")
        keys = sorted(set((k[1],k[2]) for k in tab if k[0]==d))
        for (ij,p) in keys:
            pv = Fraction(*p)
            # mark pole name
            name = {Fraction(0):'0',Fraction(1):'1',Fraction(-1):'-1',-r:'-r',r:'+r'}.get(pv,str(pv))
            val = tab[(d,ij,p)]
            print(f"  {ij} pole {name:>3}: {val}")
