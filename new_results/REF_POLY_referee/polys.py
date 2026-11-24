import mpmath as mp
def R(a, b=1): return mp.mpf(a)/b
def regular(m):
    return [(mp.cos(i*mp.pi/m), mp.sin(i*mp.pi/m)) for i in range(1, m+1)]
def affine(U, M):
    # v' = M v ; K' = M^{-T} K
    return [(M[0][0]*x+M[0][1]*y, M[1][0]*x+M[1][1]*y) for x, y in U]
def polys():
    return {
     'diamond':      regular(2),
     'parallelogram':[(R(1),R(1,2)),(R(-1,3),R(1))],
     'reg_hexagon':  regular(3),
     'affreg_hexagon': affine([(R(1),R(0)),(R(0),R(1)),(R(1),R(-1))], [[R(2),R(1,3)],[R(-1,2),R(1)]]),
     'reg_octagon':  regular(4),
     'reg_decagon':  regular(5),
     'gen_hexagon':  [(R(1),R(0)),(R(1,5),R(1)),(R(2,3),R(-1,2))],
     'gen_octagon':  [(R(1),R(0)),(R(1,5),R(1)),(R(-2,3),R(1,2)),(R(1,2),R(1,3))],
    }
