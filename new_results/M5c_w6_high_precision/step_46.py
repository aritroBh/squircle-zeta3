import pickle, json, datetime
import sys
sys.dont_write_bytecode=True
work='<project>/research/08_theory/M5c_w6_high_precision'
sys.path.insert(0,work)
import mpmath as mp
mp.mp.dps=3700
from w6_values import leaf_recipes
recs=leaf_recipes()
saved=pickle.loads(open(work+'/w6_leaves_3500.pkl','rb').read())
fix=json.load(open(work+'/mpl_recheck.json'))
key=[k for k in fix if "'mpl', 4, 2" in k]
print('key:',key)
newv=fix[key[0]][0]
idx=[i for i,r in enumerate(recs) if r[0]=='mpl' and r[1]==4 and r[2]==2][0]
print('patching index', idx, recs[idx])
old=saved[idx]
saved[idx]=newv
pickle.dump(saved, open(work+'/w6_leaves_3500.pkl','wb'))
saved2=pickle.loads(open(work+'/w6_leaves_3500.pkl','rb').read())
print('old:', old[:50])
print('new:', saved2[idx][:50])
print('diff old-new:', mp.nstr(abs(mp.mpf(old)-mp.mpf(saved2[idx])),5))
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} BUGFIX: mpl_general N-rate for |y|<1; Li_4,2(b,sqrt2/2) was good only to 1e-2582 (M5b: ~1e-1751, invalidating their 1900-digit LLL). Recomputed with rate rx*max(1,ry), radius 3.7e-3639; cache patched. Other 7 mpl leaves old-vs-new <=1.1e-3584.\n')
print('patched + logged')
