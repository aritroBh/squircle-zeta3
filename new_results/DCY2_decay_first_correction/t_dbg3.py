
from flint import ctx, arb
ctx.dps = 40
import dcy2_cert as C
N = arb(460, 60)
q = C.quantities(N)
for k in ["E_mom","E_tailmom","Tails","LHS","B3","alpha","M4","Hsup","lamL"]:
    v = q[k]
    print(k, "=", v.str(6) if hasattr(v,"str") else v)
# find where nan arises: tails
y = q['y']; Delta = q['Delta']; N0 = q['N']
from flint import arb as A
bR = C.bR; SQ2 = C.SQ2; PI = C.PI
def LogM(t):
    u = 1+y*t
    Q = u**2 + SQ2*y*u + y**2
    return (N0-1)*(u).log() + (N0/2)*(1+SQ2*y+y**2).log() - (N0/2)*(Q).log() - bR*t
LM_D = LogM(Delta); LM_mD = LogM(-Delta)
print("LM_D =", LM_D.str(6), " LM_mD =", LM_mD.str(6))
print("exp LM_D =", LM_D.exp().str(6))
uDR = 1+q['yDelta']
QR = uDR**2 + SQ2*y*uDR + y**2
LamPR = y*((N0-1)/uDR - N0*(uDR+y/SQ2)/QR) - bR
print("LamPR =", LamPR.str(6))
rDelta = min(-LamPR, bR)
print("rDelta =", rDelta.str(6))
