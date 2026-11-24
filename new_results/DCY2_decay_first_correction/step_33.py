for N in [400, 1600, 6400, 25600, 102400]:
    N = mpf(N)
    theta = sqrt(pi*N)
    T1 = -1j*sqrt(pi)/N*exp(-theta+j*theta)*(1+varpi1/sqrt(N))
    Mo = model(N)
    L = sqrt(pi)/N*exp(-theta)
    print(f"N={int(N):7d}  (model-T1)/(L/N) = {float(abs(Mo-T1)/(L/N)):10.4f}   (model-T1)/(L/sqrtN) = {float(abs(Mo-T1)/(L/sqrt(N))):.5f}")
