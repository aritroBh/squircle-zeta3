# Let's understand how words work in the engine
# The POINTS dict maps geometric points to numeric indices

# From inspection: POINTS[None] = 0, POINTS[0] = 4, POINTS[1] = 6, POINTS[-1] = 2, etc.

# Let me now check the source of find_cross and dlog
print("Source of find_cross:")
print(inspect.getsource(C.find_cross)[:2000])