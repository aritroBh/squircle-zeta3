import mpmath as mp
mp.mp.dps = 80
def P1_formula(m):
    return 4*m*mp.sin(mp.pi/(2*m))**2/mp.cos(mp.pi/(2*m))
def R1_at0(m):
    A = mp.pi/m
    V = mp.fsum([mp.sin(i*A) for i in range(1,m)])  # = cot(A/2)
    U = mp.fsum([mp.sin(i*A)*mp.log(mp.sin(i*A)) for i in range(1,m)])/V
    return mp.log(V) - U
# values from coefficient integrals
P1_int = {3: mp.mpf('2.40113226770588735429020690407746529554651722044431670825375284904281178764684126665661951'),
          4: mp.mpf('2.75023713004966371870863005296385194252854498484351236290611395868134254516571310482650575')}
print("P'(1) = P(1)*R1(0) checks:")
for m in [2,3,4]:
    pred = P1_formula(m)*R1_at0(m)
    print(f"m={m}: P(1)*R1(0) = {mp.nstr(pred,60)}")
    if m in P1_int:
        print(f"     integral value = {mp.nstr(P1_int[m],60)}")
        print(f"     residual = {mp.nstr(abs(pred-P1_int[m]),5)}")
print()
for m in [2,3,4,5,6,8]:
    print(f"m={m}: R1(0) = {mp.nstr(R1_at0(m),40)}")
# m=4 closed form check: R1(0) = log(1+sqrt2) + (1-sqrt2/2)*log2 ?
r = mp.log(1+mp.sqrt(2)) + (1-mp.sqrt(2)/2)*mp.log(2)
print("m=4 closed R1(0):", mp.nstr(r,60), " res:", mp.nstr(abs(r-R1_at0(4)),5))
