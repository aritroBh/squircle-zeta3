from cor8r_num import *
setdps(50)
# Validate d/dlog c Li_{2,2}(1/2,1-c) = Li_{2,1}(1/2,1-c) * dlog(1-c)/dlog c  numerically
def F(cval):
    return Li_ab(2,2, mp.mpf('0.5'), 1-mp.mpf(cval), terms=600)
c0 = mp.mpf('0.31')
h = mp.mpf('1e-13')
dF = (F(c0*mp.exp(h))-F(c0*mp.exp(-h)))/(2*h)
pred = Li_ab(2,1, mp.mpf('0.5'), 1-c0, terms=600) * (-c0/(1-c0))
print("dLi22/dlog c :", mp.nstr(dF,30))
print("prediction   :", mp.nstr(pred,30))
print("residual:", mp.nstr(abs(dF-pred),5))
