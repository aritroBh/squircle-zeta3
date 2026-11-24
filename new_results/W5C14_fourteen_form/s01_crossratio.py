"""Step 1 (exact): which Li5 arguments x are cross-ratios of four of the eight points
S = {0, +-1, +-a, +-b, inf} (a=1+sqrt2, b=sqrt2-1)?  Exhaustive over all ordered 4-tuples.
x is a level-8 argument iff some Mobius map phi sends (0,1,inf,x) -> four points of S,
i.e. iff x = CR(s1,s2,s3,s4) := (s4-s1)(s2-s3)/((s4-s3)(s2-s1)) (the image of s4 under the
Mobius map sending s1,s2,s3 to 0,1,inf).  Then Li5(x) = -I(0;1,0,0,0,0;x) becomes an iterated
integral between points of S with letters in S."""
import itertools, json
from qs2 import Q, A, B, R2
INF='inf'
S=[Q(0),Q(1),Q(-1),A,-A,B,-B,INF]
names=['0','1','-1','a','-a','b','-b','inf']
def mob_to_01inf(s1,s2,s3,z):
    # phi(z) = (z-s1)(s2-s3)/((z-s3)(s2-s1)), with inf handled
    def d(u,v):  # u-v, with inf -> None marker
        return None if (u is INF or v is INF) else (u-v)
    num=[d(z,s1),d(s2,s3)]; den=[d(z,s3),d(s2,s1)]
    # drop the factors containing inf (they cancel in pairs since each inf appears once in num and once in den)
    n=Q(1); m=Q(1)
    for t in num:
        if t is not None: n=n*t
    for t in den:
        if t is not None: m=m*t
    return n/m
crs={}
for t in itertools.permutations(range(8),4):
    x=mob_to_01inf(S[t[0]],S[t[1]],S[t[2]],S[t[3]])
    crs.setdefault(x,[]).append(tuple(names[i] for i in t))
print('distinct cross-ratios:',len(crs))
args={'b':B,'-b':-B,'-b^2':-(B*B),'1/sqrt2':R2/2,'1/2':Q(1,0)/2,
      '1-1/sqrt2':1-R2/2,'2-sqrt2':2-R2,
      '-b/2':-B/2,'(2-sqrt2)/4':(2-R2)/4,'(2+sqrt2)/4':(2+R2)/4,
      '2b':2*B,'(4-3sqrt2)/8':(4-3*R2)/8,'b^4':B*B*B*B,'1-b^4':1-B*B*B*B,
      'b^2':B*B,'-1':Q(-1),'-1/sqrt2':-R2/2}
out={}
for k,x in args.items():
    hit=crs.get(x)
    print(f'{k:14s} = {x!r:28s} level8={"YES" if hit else "no "}', ('e.g. '+str(hit[0])+f'  ({len(hit)} tuples)') if hit else '')
    out[k]={'x':repr(x),'level8':bool(hit),'example':hit[0] if hit else None,'ntuples':len(hit) if hit else 0}
json.dump(out,open('s01_crossratio.json','w'),indent=1)
# also: the full list of cross-ratios in (-1,1) with x and 1-x S-units (for reference)
inr=sorted([x for x in crs if abs(float(x))<1 and not x.iszero()],key=float)
print('cross-ratios in (-1,1)\\{0}:',len(inr)); print(inr)
