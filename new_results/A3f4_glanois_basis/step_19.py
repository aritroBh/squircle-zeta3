import sys, os, time, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 4200
cc.N = 1400
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 60
I = acb(0,1); T = 2*arb.pi()*I

def letter(sg): return 1 if sg==1 else 5
def zeta_word(xs, eps):
    p = len(xs); run=1; signs=[1]*p
    for k in range(p-1,-1,-1):
        run *= eps[k]; signs[k]=run
    Tt=[]
    for k in range(p):
        Tt.append(letter(signs[k])); Tt.extend([-1]*(xs[k]-1))
    return tuple(reversed(Tt))

def val(xs, eps, s=0):
    return cc.h(zeta_word(xs,eps)) * T**s

# stuffle checks (all at ~1200 digits, internal to h; validates word construction at all depths)
checks = []
# z(1;xi)z(2;xi) = z(1,2;1,xi)+z(2,1;1,xi)+z(3;xi)
checks.append(val((1,),(1,))*val((2,),(1,)) - (val((1,2),(1,1))+val((2,1),(1,1))+val((3,),(1,))))
# z(1;xi)z(2;-xi): eps product: z(1,2;1,-xi)? z(2; -1*xi8): eps=( -1,); product = z(1,2;1,-1)+z(2,1;-1,1)+z(3;-1)
checks.append(val((1,),(1,))*val((2,),(-1,)) - (val((1,2),(1,-1))+val((2,1),(-1,1))+val((3,),(-1,))))
# z(1;xi)z(1,1;1,-xi) = z(1,1,1; 1,1,-1)+z(1,1,1; 1,1,-1 perm)... full stuffle:
# z(1;a) z(1,1;b,c) = z(1,1,1; a,b,c) + z(1,1,1; b,a,c) + z(1,1,1; b,c,a) + z(1,2; b, a*c) + z(2,1; a*b, c)
checks.append(val((1,),(1,))*val((1,1),(1,-1)) - (
    val((1,1,1),(1,1,-1)) + val((1,1,1),(1,1,-1))*0 + val((1,1,1),(1,1,-1))*0 ))
# do the general one with explicit terms:
a=(1,); bc=(1,-1)
full = ( val((1,1,1),(1,1,-1))   # a before b,c
       + val((1,1,1),(1,1,-1)) ) # placeholder
# write it properly: z(1;eps_a) * z(1,1;eps_b,eps_c):
#  = z(1,1,1; ea,eb,ec) + z(1,1,1; eb,ea,ec) + z(1,1,1; eb,ec,ea)
#    + z(2,1; ea*eb, ec) + z(1,2; eb, ea*ec)
ea,eb,ec = 1, 1, -1
lhs = val((1,),(ea,))*val((1,1),(eb,ec))
rhs = ( val((1,1,1),(ea,eb,ec)) + val((1,1,1),(eb,ea,ec)) + val((1,1,1),(eb,ec,ea))
      + val((2,1),(ea*eb,ec)) + val((1,2),(eb,ea*ec)) )
checks.append(lhs-rhs)
worst=mp.mpf(0)
for d in checks:
    e = abs(cc.approx(d)); worst=max(worst,e)
    print('residual:', mp.nstr(e,3))
print('worst:', mp.nstr(worst,3))
