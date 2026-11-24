from towlib import *
def nu(i,j,lvl,delta):
    return u(i,j,lvl,delta) - u(i,j,lvl-1,delta)

for key in [(2,1),(3,0),(2,2)]:
    i,j = key; d = i+j
    print(f'--- key {key}, d={d} ---')
    for delta in range(0, 2*d-3):
        row = []
        for k in range(d, 14):
            row.append(u(i,j,k,delta))
        print(f'  delta={delta}: u(k) for k={d}..13:')
        print('   ', [str(x) for x in row])
