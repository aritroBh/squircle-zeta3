"""AD identities for delta<=2 with TOW2 machinery (TOW1's E1,E2a,E2b re-derived in TOW2 conventions)."""
import sys; sys.path.insert(0,'.')
import prove_adj34 as A
for (delta,s) in [(0,0),(1,0),(2,0),(2,1)]:
    A.main(delta, s)
