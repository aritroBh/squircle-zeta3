from relations import *
from check_conversion import evaluate, approx, mp

def main():
    for n in range(1,5):
        bad=[]; total=0; worst=mp.mpf(0)
        if n==1:
            rows,counts=make_relations(n,{})
            candidates=[('seed',dict(r)) for r in rows]
        else:
            candidates=[('DS',(ds(u,v))) for k in range(1,n) for u in words(k,False) for v in words(n-k)]
        candidates += [(f'D{d}',r) for d in (2,4,8) for r in distribution(n,d)]
        for name,row in candidates:
            if not row: continue
            err=abs(approx(evaluate(row))); total+=1; worst=max(worst,err)
            if err>mp.mpf('1e-30'):
                bad.append((name,row,mp.nstr(err,10)))
                if len(bad)<=5: print('BAD',n,bad[-1],flush=True)
        log(f'Relation numeric checks weight={n}: total={total}, bad={len(bad)}, worst={mp.nstr(worst,10)}')
        assert not bad

if __name__=='__main__': main()
