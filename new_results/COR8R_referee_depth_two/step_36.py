from cor8r_num import *
setdps(40)
c0 = mp.mpf('0.31')
A = 1/(1+c0); B = 2/(1+c0)
def Gv(word, z=mp.mpf(1)):
    # numeric regularized G(word;1): repeated/trailing zeros handled via shuffles at value level
    return G_int(list(word), z)
# check: shuffle((B,),(A,A)): G(B;1)*G(A,A;1) = G(B,A,A)+G(A,B,A)+G(A,A,B)
GB = mp.log(1-1/B)
GA = mp.log(1-1/A)
GAA = GA*GA/2
GBAA = G_int([B,A,A],1)
GABA = G_int([A,B,A],1)
GAAB = G_int([A,A,B],1)
lhs = GB*GAA
rhs = GBAA+GABA+GAAB
print("shuffle (B)x(A,A): lhs=", mp.nstr(lhs,25), " rhs=", mp.nstr(rhs,25), " diff:", mp.nstr(abs(lhs-rhs),5))
# now the bad one: u=(A,), v=(B,0): G(A;1)G(B,0;1) = G(A,B,0)+G(B,A,0)+G(B,0,A)
# G(B,0;1) = -G(0,B;1) = Li_2(1/B);  G(B,0;1) numeric: int dt/(t-B) log t
GB0 = G_int([B,0],1)  # should work: inner log t integrable
print("G(B,0;1) numeric:", mp.nstr(GB0,25), " Li2(1/B):", mp.nstr(Li_n(2,1/B),25))
GBA0 = G_int([B,A,0],1)
GB0A = G_int([B,0,A],1)
GAB0 = G_int([A,B,0],1)
lhs = GA*GB0
rhs = GAB0+GBA0+GB0A
print("shuffle (A)x(B,0): lhs=", mp.nstr(lhs,25), " rhs=", mp.nstr(rhs,25), " diff:", mp.nstr(abs(lhs-rhs),5))
