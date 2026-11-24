mp.mp.dps = 80

# Tests
print("Testing surface area calculation:")
S2 = surface_area_correct(2)
print(f"S(2) = {S2}")
print(f"Expected 4*pi = {4*mp.pi}")
print(f"Relative error: {abs(S2 - 4*mp.pi)/(4*mp.pi)}")
print()

# Test p = 1 (octahedron): S(1) should be 4*sqrt(3)
# Note: our function may have issues at p = 1, let's test near p = 1
for p_test in [1.001, 1.01, 1.1]:
    S_test = surface_area_correct(p_test)
    print(f"S({p_test}) = {S_test}")

print("\nExpected S(1) = 4*sqrt(3) =", 4*mp.sqrt(3))