
from flint import arb, ctx
ctx.dps = 30
for m,r in [(1.205,0.1224),(1.0826,0.1224),(10.0,0.5),(100.0,1.0),(0.0,0.1224),(1.205,0.0),(5.5,2.25),(1e6,1e4)]:
    e = arb(m,r)
    print(f"arb({m},{r}): str={e}  mid={e.mid()} rad={e.rad()}")
