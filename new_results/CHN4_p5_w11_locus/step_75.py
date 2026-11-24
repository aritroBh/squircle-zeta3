import pickle, numpy as np
d=pickle.load(open('signmap.pkl','rb')); grid=d['grid']; signmap=d['signmap']
N=len(grid)
# ASCII map: rows=a1 (top=large), cols=a2. '+'=only pos, '-'=only neg, 'B'=both, '.'=no pts
print("Sign map of Q11 on V10+ (a1 vertical, a2 horizontal, log [1e-3,1e3]):")
print("    a2 ->")
hdr="     "+"".join(str((j//10)%10) for j in range(0,N,3))
print(hdr)
for i in range(N-1,-1,-3):
    row=""
    for j in range(0,N,3):
        s=signmap.get((i,j))
        if s is None: row+='.'
        elif s=={1}: row+='+'
        elif s=={-1}: row+='-'
        else: row+='B'
    lbl=f"a1={grid[i]:7.2g}" if i%6==0 else "        "
    print(f"{lbl} {row}")
print("\nLegend: +=Q11>0 only, -=Q11<0 only, B=both signs (Gamma crosses), .=no V10+ pt")
# Count structure
both=[(i,j) for (i,j),s in signmap.items() if len(s)>1]
print(f"\nBoth-sign cells: {len(both)} spread i in [{min(b[0] for b in both)},{max(b[0] for b in both)}], j in [{min(b[1] for b in both)},{max(b[1] for b in both)}]")
