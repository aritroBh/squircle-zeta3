"""Standalone exact verifier for the M4b W4'' certificate (adapted from A3h5 verify_D.py).

Rebuilds EVERY row from first principles; reads none of the M4b pipeline pickles:
  * T0 rows: relations.make_relations regenerated; lower weight bases asserted subset of generated;
  * frame rows: re-derived from A3h4 provenance via candidates/convert/identity, asserted equal to exports;
  * atom rows: classical weight-4 constants rebuilt from the certificate's atom list via cmzv.atom/mul;
  * target W4'': rebuilt from M4a_decomposition/gpl_w4pp.tsv (parsed here), cross-checked against the
    omega4.json pole data (gexpr logic), converted via cmzv dlog/forms_path exactly as cmzv.vectors().
Then checks, with exact Fractions on all 5184 admissible weight-4 level-8 words:
    W4''_X = sum_k c_k^X B_k + sum_j z_j^X F_j + sum_i y_i^X T0_i ,  X in {rat, sq}.
Hence W4'' = R'' = sum_k (c_k^rat + sqrt2 c_k^sq) B_k as real numbers (T0 and frame rows vanish by
regularised double shuffle + distribution and by function-level double shuffle).
"""
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/A3h4_frame_ds')  # ADAPTED
sys.path.insert(0, '<project>/research/08_theory/A3c_claude_elementary/Q')  # ADAPTED

from frame_ds import *

import time as _time  # ADAPTED: log into the M4b folder
def log(s):
    s=_time.strftime('%Y-%m-%d %H:%M:%S ')+s
    print(s, flush=True)
    with open('<project>/research/08_theory/M4b_w4pp_exact/verifier.log','a') as _f: _f.write(s+'\n')

import re as _re
import cmzv as cm
import relations as standard

M4B = '<project>/research/08_theory/M4b_w4pp_exact'  # ADAPTED
M4A = '<project>/research/08_theory/M4a_decomposition'  # ADAPTED

def rebuild_standard():
    standard.log=log  # All logs stay in this work directory.
    selected={}; maps={}
    for n in (1,2,3):
        generated,_=standard.make_relations(n,selected)
        oldcols,oldrows,basis=load(THEORY/f'A3c_claude_elementary/Q/weight{n}.pickle')
        assert tuple(oldcols)==words(n)
        generated=set(generated)
        assert all(row in generated for row in basis),('lower basis not regenerated',n)
        selected[n]=basis
        maps[n]=basis
    generated,_=standard.make_relations(4,{})
    return generated,maps

def parse_qs2(s):
    s=s.strip()
    m=_re.fullmatch(r"\((.+?)([+-])(\d+(?:/\d+)?)\*sqrt2\)",s)
    if m:
        p=F(m.group(1)); q=F(m.group(3))
        if m.group(2)=='-': q=-q
        return Q(p,q)
    m=_re.fullmatch(r"(-?\d+(?:/\d+)?)\*sqrt2",s)
    if m: return Q(0,F(m.group(1)))
    return Q(F(s))

def parse_word(s):
    s=s.strip()
    assert s.startswith('(') and s.endswith(')')
    body=s[1:-1].strip()
    if not body: return ()
    parts,depth,cur=[],0,''
    for ch in body:
        if ch=='(': depth+=1; cur+=ch
        elif ch==')': depth-=1; cur+=ch
        elif ch==',' and depth==0: parts.append(cur); cur=''
        else: cur+=ch
    parts.append(cur)
    return tuple(parse_qs2(x) for x in parts if x.strip())

def tsv_terms():
    rows=[]
    with open(M4A+'/gpl_w4pp.tsv') as f:
        for line in f:
            if line.startswith('#') or line.startswith('coefficient_p'): continue
            p,q,k,w=line.rstrip('\n').split('\t')
            rows.append((parse_word(w),int(k),Q(F(p),F(q))))
    out={}
    for w,k,c in rows: out[(w,k)]=out.get((w,k),Q(0))+c
    return out

def gexpr_terms():
    """Independent re-derivation from omega4.json (gexpr logic; exact IBP pole data)."""
    import w4_symbol as ws
    om=json.load(open(M4A+'/omega4.json'))
    om={tuple(map(int,kk.split(','))):v for kk,v in om.items()}
    POLE={'0':Q(0),'a':cm.A,'-a':-cm.A,'b':cm.B,'-b':-cm.B}
    out={}
    for (i,j) in ((2,1),(3,0)):
        for name,(p,q) in om[(i,j)]:
            z=POLE[name]; coefficient=Q(F(p),F(q))
            for (w,k),d in ws.log_product(*((ws.LT,)*i+(ws.MT,)*j)).items():
                out[((z,)+w,k)]=out.get(((z,)+w,k),Q(0))+coefficient*d
    return {kk:v for kk,v in out.items() if v}

def rebuild_target():
    tsv=tsv_terms()
    assert tsv==gexpr_terms(),'gpl_w4pp.tsv does not match omega4.json re-derivation'
    log(f'target: gpl_w4pp.tsv == omega4.json re-derivation ({len(tsv)} terms); weight check '
        f'{all(len(w)+k==4 for (w,k) in tsv)}')
    lhs=[{},{}]
    for (w,k),coefficient in tsv.items():
        forms=[dlog(p) for p in w]
        term=add((forms_path(forms,POINTS[Q(1)],POINTS[cm.A]),1),
                 (forms_path(forms,POINTS[Q(1)],POINTS[cm.B]),-1))
        for _ in range(k): term=mul(term,cm.atom(('log',Q(2))))
        for j,c in enumerate((coefficient.p,coefficient.q)):
            lhs[j]=add((lhs[j],1),(term,c))
    return lhs

def rebuild_atom(enc):
    ats=[]
    for a in enc:
        if a[0]=='log': ats.append(('log',parse_qs2(a[1])))
        elif a[0] in ('pi2','z3'): ats.append((a[0],))
        elif a[0]=='li': ats.append(('li',int(a[1]),parse_qs2(a[2])))
        else: raise ValueError(a)
    v={():F(1)}
    for a in ats: v=mul(v,cm.atom(a))
    return v

def main():
    cert=json.loads(open(M4B+'/W4pp_certificate.json').read())  # ADAPTED
    std,bases=rebuild_standard()
    rebuilt={}; counts=Counter()
    for label in cert['labels']:
        key=repr(label)
        if label[0]=='T0':
            lab=label[1]; family=lab[0]
            if family=='std': row=std[lab[1]]
            else:
                assert family in ('lift1','lift2old','lift3old'),lab
                n={'lift1':1,'lift2old':2,'lift3old':3}[family]
                row=primitive(mul(dict(bases[n][lab[1]]),{tuple(lab[2]):1}))
            counts[family]+=1
        elif label[0]=='frame':
            _,n,i,w=label; data=load(HERE/f'validated_w{n}.pickle')
            g,fs=data['provenance'][i]
            p0,p1,pinf=g
            ts={cross(p,p0,p1,pinf):p for p in POINTS}
            bs={1/t for t in ts if t is not None and abs(t)>1}
            # Regenerate and validate the complete stuffle grid in exact Q(sqrt2).
            assert fs in tuple(candidates(n,bs,len(fs)))
            raw=primitive(convert(g,ts,identity(fs)))
            assert raw==data['rows'][i],('frame export mismatch',n,i)
            row=primitive(mul(dict(raw),{tuple(w):1}))
            counts[f'frame_w{n}']+=1
        else:
            assert label[0]=='atom'
            row=tuple(sorted(rebuild_atom(cert['atoms'][label[1]]).items()))
            counts['atom']+=1
        rebuilt[key]=dict(row)
    targets=rebuild_target()
    ok=[]
    for j,(target,coeffs) in enumerate(zip(targets,cert['coefficients'])):
        residual=Counter({w:-a for w,a in target.items()})
        for label,a in zip(cert['labels'],coeffs):
            a=F(a)
            for w,b in rebuilt[repr(label)].items(): residual[w]+=a*b
        nonzero={w:a for w,a in residual.items() if a}
        log(f'EXACT rebuilt W4pp component={j} support={len(target)} residual_nonzero={len(nonzero)}')
        assert not nonzero,('exact identity failed',j,next(iter(nonzero.items()),None))
        ok.append(not nonzero)
    log(f'EXACT certificate source families={dict(counts)} total_used_rows={len(cert["labels"])}')
    # readable R''
    names=cert['atom_names']; cr=[F(a) for a in cert['coefficients'][0]]; cs=[F(a) for a in cert['coefficients'][1]]
    used=[(names[l[1]],cr[i],cs[i]) for i,l in enumerate(cert['labels']) if l[0]=='atom']
    log("R'' = sum over classical atoms of (c_rat + c_sq*sqrt2)*atom:")
    for n,a,b in used: log(f'   ({a} {"+" if b>=0 else "-"} {abs(b)}*sqrt2) * {n}')
    log('VERDICT: ' + ("W4'' = R'' PROVED" if all(ok) else 'NOT established'))

if __name__=='__main__': main()
