import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
src=open('mzv_exact.py').read()
old="""def rref_flint(rels, idx):
    from flint import fmpq_mat, fmpq
    n=len(idx)
    rows=[]
    for r in rels:
        rows.append((r,))
    M=fmpq_mat(len(rels), n, [c for (r,) in rows for c in [r.get(w,0) for w in idx]])
    R, rank = M.rref()"""
new="""def rref_flint(rels, idx):
    from flint import fmpq_mat, fmpq
    n=len(idx); wlist=list(idx)
    data=[]
    for r in rels:
        for w in wlist:
            c=r.get(w)
            data.append(fmpq(c.numerator,c.denominator) if c is not None else fmpq(0))
    M=fmpq_mat(len(rels), n, data)
    R, rank = M.rref()"""
src=src.replace(old,new)
open('mzv_exact.py','w').write(src)
t0=time.time()
import mzv_exact as M
T10=M.DSTable(10)
T10f=M.DSTableFlint(10)
print('nullities:', T10.nullity, T10f.nullity)
same=all(T10.reduce_word(w)==T10f.reduce_word(w) for w in M.WORDS[10])
print('W=10 reductions identical:', same)
t=time.time()
T12=M.DSTableFlint(12)
print('W=12 nullity:', T12.nullity, 'expected 12', '(%.0fs)'%(time.time()-t))
# validations at W=12
dualok=all(T12.reduce_word(w)==T12.reduce_word(M.dual(w)) for w in M.WORDS[12])
print('duality all weight-12 words:', dualok)
s2={}
for a in range(2,12): s2[(a,12-a)]=F(1)
print('sum theorem depth2 W12:', T12.reduce(s2)==[F(1)]+[F(0)]*11)
print('z(2,1^10)=z(12):', T12.reduce_word((2,)+(1,)*10)==[F(1)]+[F(0)]*11)
print('t=%.1fs'%(time.time()-t0))
