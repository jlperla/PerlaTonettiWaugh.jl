#=
    Rescaled discretization tests. 
=#

ξ_1 = 1.0 
ξ_2 = 2.0
x_uniform = 1:5 
x_irregular = [1, 2, 3, 4, 5] # This is an AbstractArray and not a Range, so it calls the right method. 

# Test for uniform grid code, ξ_1.

σ = 1; μ = -1;
x, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(x_uniform,ξ_1) # Dispatches on the discrete code. 
row1 = [-1 + (1+ξ_1) + (-2+1+ξ_1)/2 , 0.5, 0.0, 0.0, 0.0]'
row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_1)/2]'
@test μ * L_1_minus + σ^2/2 * L_2 == cat(1, row1, row2, row3, row4, row5) # only test for backwards now
# Test for irregular grid code, ξ_1. 
# Test irregular grid function produce proper regular grid
σ = 1; μ = -1;
x, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(x_irregular,ξ_1) # Dispatches on the discrete code. 
row1 = [-1 + (1+ξ_1) + (-2+1+ξ_1)/2 , 0.5, 0.0, 0.0, 0.0]'
row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_1)/2]'
@test μ * L_1_minus + σ^2/2 * L_2 == cat(1, row1, row2, row3, row4, row5) # only test for backwards now
# Uniform ... ξ_2. 
σ = 1; μ = -1;
x, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(x_uniform,ξ_2) # Dispatches on the discrete code. 
row1 = [-1 + (1+ξ_2) + (-2+1+ξ_2)/2 , 0.5, 0.0, 0.0, 0.0]'
row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_2)/2]'
@test μ * L_1_minus + σ^2/2 * L_2 == cat(1, row1, row2, row3, row4, row5) # only test for backwards now
# Irregular ... ξ_2.  
σ = 1; μ = -1;
x, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(x_irregular,ξ_2) # Dispatches on the discrete code. 
row1 = [-1 + (1+ξ_2) + (-2+1+ξ_2)/2 , 0.5, 0.0, 0.0, 0.0]'
row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_2)/2]'
@test μ * L_1_minus + σ^2/2 * L_2 == cat(1, row1, row2, row3, row4, row5) # only test for backwards now

# Test for error handling on ξ (i.e., if ξ can't be negative or is bounded or something like that) for each method. 
# Test for proper dispatch (i.e., type of the first return from uniform is a Range).


#=
    ODE Rescaling tests.
=#

    # Parameters
    const ξ = 1.0 # Constant rescaling value. 
    const r = 0.05 - 0.021 # Constant state multiplier, this is picked close to r-g(t) in the growth model(residualtest.jl)
    const μ_val2 = 0.005 - 0.021 # Constant drift < 0. This is picked close to γ-g(t) in the growth model
    const σ_val2 = 0.02 # Constant shock, same, picked close to growth model
    const π(ξ) = z -> exp(ξ * z) # Rescaled π function. This defines a family of π functions, parametrized by ξ. 

    # Grid 
    x_min = 0.0
    x_max = 5.0
    M = 200
    x = linspace(x_min,x_max,M)# Some partition of [0, 5].

    #Tests for stability of irregular grid
    x_irregular=unique([linspace(x_min, 1.0, 500)' linspace(1.0, x_max, 201)'])
    # add point in the bottom
    x_irregular_add1=unique([linspace(x_min, 1.0, 700)' linspace(1.0, x_max, 201)'])
    # add point in the top
    x_irregular_add2=unique([linspace(x_min, 1.0, 500)' linspace(1.0, x_max, 301)'])
    # change x_max
    x_irregular_add3=unique([linspace(x_min, 1.0, 500)' linspace(1.0, x_max+1.0, 201)'])

    #Rescaling method
    x_irregular, L_1_tild, L_1_plus_tild, L_2_tild = rescaled_diffusionoperators(x_irregular, ξ)
    # modify r and μ_val2 ,π to match (13)
    r_tild = r - ξ*μ_val2 - σ_val2^2/2*ξ^2;
    μ_tild = μ_val2 + σ_val2^2*ξ;
    π_tild(x) = 1;
    L_T_tild_irregular = r_tild*I - μ_tild*L_1_tild - σ_val2^2/2 * L_2_tild # Construct the aggregate operator. 
    v_tild_irregular = L_T_tild_irregular \ π_tild.(x_irregular) # Solution to the rescaled differential equation. 
    # Use addpoint Grid
    x_irregular_add1, L_1_tild, L_1_plus_tild, L_2_tild = rescaled_diffusionoperators(x_irregular_add1, ξ)
    L_T_tild_irregular_add1 = r_tild*I - μ_tild*L_1_tild - σ_val2^2/2 * L_2_tild # Construct the aggregate operator. 
    v_tild_irregular_add1 = L_T_tild_irregular_add1 \ π_tild.(x_irregular_add1)
    # Use addpoint grid in top
    x_irregular_add2, L_1_tild, L_1_plus_tild, L_2_tild = rescaled_diffusionoperators(x_irregular_add2, ξ)
    L_T_tild_irregular_add2 = r_tild*I - μ_tild*L_1_tild - σ_val2^2/2 * L_2_tild # Construct the aggregate operator. 
    v_tild_irregular_add2 = L_T_tild_irregular_add2 \ π_tild.(x_irregular_add2)
    # Use change x_max grid 
    x_irregular_add3, L_1_tild, L_1_plus_tild, L_2_tild = rescaled_diffusionoperators(x_irregular_add3, ξ)
    L_T_tild_irregular_add3 = r_tild*I - μ_tild*L_1_tild - σ_val2^2/2 * L_2_tild # Construct the aggregate operator. 
    v_tild_irregular_add3 = L_T_tild_irregular_add3 \ π_tild.(x_irregular_add3)

    # Interpolaration
    v_int1=LinInterp(x_irregular_add1, v_tild_irregular_add1)
    v_int2=LinInterp(x_irregular_add2, v_tild_irregular_add2)
    v_int3=LinInterp(x_irregular_add3, v_tild_irregular_add3)

    @test norm(v_int1.(x_irregular)-v_tild_irregular,Inf)≈ 0.0 atol = 1e-2
    @test norm(v_int2.(x_irregular)-v_tild_irregular,Inf)≈ 0.0 atol = 1e-2
    @test norm(v_int3.(x_irregular)-v_tild_irregular,Inf)≈ 0.0 atol = 3e-1

    plot(x_irregular,v_tild_irregular)
    plot!(x_irregular_add1,v_tild_irregular_add1)
    plot!(x_irregular_add2,v_tild_irregular_add2)
    plot!(x_irregular_add3,v_tild_irregular_add3)

    v_rescale_irr = v_tild_irregular.*exp.(ξ*x_irregular)
    v_rescale_irr_add1 = v_int1.(x_irregular).*exp.(ξ*x_irregular)
    v_rescale_irr_add2 = v_int2.(x_irregular).*exp.(ξ*x_irregular)
    v_rescale_irr_add3 = v_int3.(x_irregular).*exp.(ξ*x_irregular)

    @test norm(v_rescale_irr-v_rescale_irr_add1,Inf)≈ 0.0 atol = 1e-2

    # For the 3rd case
    indx=maximum(find(x_irregular.<=4))
    @test norm(v_int2.(x_irregular[1:indx])-v_tild_irregular[1:indx],Inf)≈ 0.0 atol = 1e-2
    @test norm(v_int3.(x_irregular[1:indx])-v_tild_irregular[1:indx],Inf)≈ 0.0 atol = 1e-2
    
