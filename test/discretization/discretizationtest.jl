using PerlaTonettiWaugh, Base.Test, BandedMatrices

include("../../src/irregulardiffusionoperators.jl")
include("../../src/diffusionoperators.jl")

## REGULAR DISCRETIZATION


## IRREGULAR DISCRETIZATION
# Test that we properly generalize the regular case for negative drift 
σ = 1; μ = -1; x = [1, 2, 3, 4, 5]
x, L_1_minus, L_1_plus, L_2 = irregulardiffusionoperators(x, 5)
row1 = [-0.5, 0.5, 0.0, 0.0, 0.0]'
row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
row5 = [0.0, 0.0, 0.0, 1.5, -1.5]'
@test μ * L_1_minus + σ^2/2 * L_2 == cat(1, row1, row2, row3, row4, row5)

# Test that we properly mirror the MATLAB on irregular grids and negative drift
x = [1.1, 2.2, 3.5, 3.7, 4.0]
x, L_1_minus, L_1_plus, L_2 = irregulardiffusionoperators(x, 5)
row1 = [-0.4132    0.4132         0         0         0]
row2 = [1.2879   -1.6084    0.3205         0         0]
row3 = [ 0    1.2821   -4.6154    3.3333         0]
row4 = [0         0   15.0000  -21.6667    6.6667]
row5 = [0         0         0    8.8889   -8.8889]
@test μ * L_1_minus + σ^2/2 * L_2 ≈ cat(1, row1, row2, row3, row4, row5) atol = 1e-4

# Test that this works even for very finely spaced grids 

# Test that this works for....

