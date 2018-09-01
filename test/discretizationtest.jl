@testset "ODE scaling" begin
#= 
Test on reascaling method and make comparison
=#
    # Parameters
     ξ = 1.0 # Constant rescaling value. 
     r = 0.05 - 0.021 # Constant state multiplier, this is picked close to r-g(t) in the growth model(residualtest.jl)
     μ_val2 = 0.005 - 0.021 # Constant drift < 0. This is picked close to γ-g(t) in the growth model
     σ_val2 = 0.02 # Constant shock, same, picked close to growth model
     π(ξ) = z -> exp(ξ * z) # Rescaled π function. This defines a family of π functions, parametrized by ξ. 

# Grid 
    z_min = 0.0
    z_max = 5.0
    M = 200
    z = range(z_min, stop = z_max, length = M)# Some partition of [0, 5].

# Vanilla (old) solution 
    z,  L_1, L_1_plus, L_2 = diffusionoperators(z) # Regular code. 
    L_T = r*I - μ_val2*L_1 - σ_val2^2/2 * L_2 # Construct the aggregate operator. 
    v = L_T \ π(ξ).(z) # Solution to the differential equation. 

# New solution. 
    z,  L_1_tilde, L_1_plus_tilde, L_2_tilde = rescaled_diffusionoperators(z,  ξ)
    # modify r and μ_val2 ,π to match (13)
    r_tilde = r - ξ*μ_val2 - (σ_val2^2/2)*ξ^2;
    μ_tilde = μ_val2 + σ_val2^2*ξ;
    π_tilde = x -> 1;
    L_T_tilde = r_tilde*I - μ_tilde*L_1_tilde - σ_val2^2/2 * L_2_tilde # Construct the aggregate operator. 
    v_tilde = L_T_tilde \ π_tilde.(z) # Solution to the rescaled differential equation.        

# Tests for uniform grid
    # Do the conversion
    v_rescale = v_tilde .* exp.(ξ*z)
    # Check for absolute similarity (norms, etc.)
    @test_broken norm(v-v_rescale,Inf) ≈ 0.0 atol = 1e-2

    # Check for greater stability with the rescaled problem. 
    z_alt = range(z_min, stop = z_max, length = M+100)
    z_alt, L_1_tilde, L_1_plus_tilde, L_2_tilde = rescaled_diffusionoperators(z_alt, ξ)
    # modify r and μ_val2 ,π to match (13)
    r_tilde = r - ξ*μ_val2 - σ_val2^2/2*ξ^2;
    μ_tilde = μ_val2 + σ_val2^2*ξ;
    π_tilde(z) = 1;
    L_T_tilde = r_tilde*I - μ_tilde*L_1_tilde - σ_val2^2/2 * L_2_tilde # Construct the aggregate operator. 
    v_tilde_alt = L_T_tilde \ π_tilde.(z_alt) # Solution to the rescaled differential equation. 
    
    # test the pre-conversion difference of v_tilde
    v_int_tilde=LinInterp(z,  v_tilde)
    @test norm(v_tilde_alt-v_int_tilde.(z_alt)) ≈ 0.0 atol = 3e-1

    # conversion to v
    v_rescale_alt = v_tilde_alt.*exp.(ξ*z_alt)
    v_int=LinInterp(z,  v_rescale)
    @test_broken norm(v_rescale_alt-v_int.(z_alt), Inf) ≈ 0.0 atol = 2e-1

    # use the ordinary method
    z_alt, L_1, L_1_plus, L_2 = diffusionoperators(z_alt) # Regular code. 
    L_T = r*I - μ_val2*L_1 - σ_val2^2/2 * L_2 # Construct the aggregate operator. 
    v_alt = L_T \ π(ξ).(z_alt) # Solution to the differential equation. 
    v_int_old=LinInterp(z,  v)
    @test_broken norm(v_alt-v_int_old.(z_alt), Inf) ≈ 0.0 atol = 2e-1

#Tests for stability of irregular grid
    z_irregular=unique([range(z_min, stop = 1.0, length = 500)' range(1.0, stop = z_max, length = 201)'])
    # add point in the bottom
    z_irregular_add1=unique([range(z_min, stop = 1.0, length = 700)' range(1.0, stop = z_max, length = 201)'])
    # add point in the top
    z_irregular_add2=unique([range(z_min, stop = 1.0, length = 500)' range(1.0, stop = z_max, length = 301)'])
    # change z_max
    z_irregular_add3=unique([range(z_min, stop = 1.0, length = 500)' range(1.0, stop = z_max+1.0, length = 201)'])

    #Rescaling method
    z_irregular, L_1_tilde, L_1_plus_tilde, L_2_tilde = rescaled_diffusionoperators(z_irregular, ξ)
    # modify r and μ_val2 ,π to match (13)
    r_tilde = r - ξ*μ_val2 - σ_val2^2/2*ξ^2;
    μ_tilde = μ_val2 + σ_val2^2*ξ;
    π_tilde(z) = 1;
    L_T_tilde_irregular = r_tilde*I - μ_tilde*L_1_tilde - σ_val2^2/2 * L_2_tilde # Construct the aggregate operator. 
    v_tilde_irregular = L_T_tilde_irregular \ π_tilde.(z_irregular) # Solution to the rescaled differential equation. 
    # Use addpoint Grid
    z_irregular_add1, L_1_tilde, L_1_plus_tilde, L_2_tilde = rescaled_diffusionoperators(z_irregular_add1, ξ)
    L_T_tilde_irregular_add1 = r_tilde*I - μ_tilde*L_1_tilde - σ_val2^2/2 * L_2_tilde # Construct the aggregate operator. 
    v_tilde_irregular_add1 = L_T_tilde_irregular_add1 \ π_tilde.(z_irregular_add1)
    # Use addpoint grid in top
    z_irregular_add2, L_1_tilde, L_1_plus_tilde, L_2_tilde = rescaled_diffusionoperators(z_irregular_add2, ξ)
    L_T_tilde_irregular_add2 = r_tilde*I - μ_tilde*L_1_tilde - σ_val2^2/2 * L_2_tilde # Construct the aggregate operator. 
    v_tilde_irregular_add2 = L_T_tilde_irregular_add2 \ π_tilde.(z_irregular_add2)
    # Use change z_max grid 
    z_irregular_add3, L_1_tilde, L_1_plus_tilde, L_2_tilde = rescaled_diffusionoperators(z_irregular_add3, ξ)
    L_T_tilde_irregular_add3 = r_tilde*I - μ_tilde*L_1_tilde - σ_val2^2/2 * L_2_tilde # Construct the aggregate operator. 
    v_tilde_irregular_add3 = L_T_tilde_irregular_add3 \ π_tilde.(z_irregular_add3)

    # Interpolation
    v_int1=LinInterp(z_irregular_add1, v_tilde_irregular_add1)
    v_int2=LinInterp(z_irregular_add2, v_tilde_irregular_add2)
    v_int3=LinInterp(z_irregular_add3, v_tilde_irregular_add3)

    @test norm(v_int1.(z_irregular)-v_tilde_irregular,Inf)≈ 0.0 atol = 1e-2
    @test norm(v_int2.(z_irregular)-v_tilde_irregular,Inf)≈ 0.0 atol = 1e-2
    @test norm(v_int3.(z_irregular)-v_tilde_irregular,Inf)≈ 0.0 atol = 3e-1

    v_rescale_irr = v_tilde_irregular.*exp.(ξ*z_irregular)
    v_rescale_irr_add1 = v_int1.(z_irregular).*exp.(ξ*z_irregular)
    v_rescale_irr_add2 = v_int2.(z_irregular).*exp.(ξ*z_irregular)
    v_rescale_irr_add3 = v_int3.(z_irregular).*exp.(ξ*z_irregular)

    @test norm(v_rescale_irr-v_rescale_irr_add1,Inf)≈ 0.0 atol = 1e-2

    # For the 3rd case
    indx=maximum(findall(z_irregular.<=4))
    @test norm(v_int2.(z_irregular[1:indx])-v_tilde_irregular[1:indx],Inf)≈ 0.0 atol = 1e-2
    @test norm(v_int3.(z_irregular[1:indx])-v_tilde_irregular[1:indx],Inf)≈ 0.0 atol = 1e-2
end 

@testset "Operator discretization" begin
# Rescaled
    ξ_1 = 1.0 
    ξ_2 = 2.0
    z_uniform = 1:5 
    z_irregular = [1, 2, 3, 4, 5] # This is an AbstractArray and not a AbstractRange, so it calls the right method. 

    # Test for uniform grid code, ξ_1.
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z_uniform,ξ_1) # Dispatches on the discrete code. 
    row1 = [-1 + (1+ξ_1) + (-2+1+ξ_1)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_1)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1) # only test for backwards now
    # Test for irregular grid code, ξ_1. 
    # Test irregular grid function produce proper regular grid
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z_irregular,ξ_1) # Dispatches on the discrete code. 
    row1 = [-1 + (1+ξ_1) + (-2+1+ξ_1)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_1)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1) # only test for backwards now
    # Uniform ... ξ_2. 
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z_uniform,ξ_2) # Dispatches on the discrete code. 
    row1 = [-1 + (1+ξ_2) + (-2+1+ξ_2)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_2)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1) # only test for backwards now
    # Irregular ... ξ_2.  
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z_irregular,ξ_2) # Dispatches on the discrete code. 
    row1 = [-1 + (1+ξ_2) + (-2+1+ξ_2)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_2)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1) # only test for backwards now

    # Test for error handling on ξ (i.e., if ξ can't be negative or is bounded or something like that) for each method. 
    # Test for proper dispatch (i.e., type of the first return from uniform is a AbstractRange).
    
# Unscaled. 
    ## REGULAR DISCRETIZATION
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = diffusionoperators(1:5) # Dispatches on the discrete code. 
    row1 = [-0.5, 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1.5]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1)
    @test -μ * L_1_plus + σ^2/2 * L_2 == [-1.5 1.5 0.0 0.0 0.0; 0.5 -2.0 1.5 0.0 0.0; 0.0 0.50 -2.0 1.50 0.0; 0.0 0.0 0.50 -2.0 1.50; 0.0 0.0 0.0 0.50 -0.50] # Test for positive drift. 

    ## IRREGULAR DISCRETIZATION
    # Test that we properly generalize the regular case 
    z = [1, 2, 3, 4, 5] # This is an AbstractArray and not a AbstractRange, so it calls the right method. 
    z,  L_1_minus, L_1_plus, L_2 = diffusionoperators(z)
    row1 = [-0.5, 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1.5]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1)
    @test -μ * L_1_plus + σ^2/2 * L_2 == [-1.5 1.5 0.0 0.0 0.0; 0.5 -2.0 1.5 0.0 0.0; 0.0 0.50 -2.0 1.50 0.0; 0.0 0.0 0.50 -2.0 1.50; 0.0 0.0 0.0 0.50 -0.50] # Test for positive drift. 

    # Test that we properly mirror the MATLAB on strictly irregular grids 
    z = [1.1, 2.2, 3.5, 3.7, 4.0]
    z,  L_1_minus, L_1_plus, L_2 = diffusionoperators(z)
    row1 = [-0.4132    0.4132         0         0         0]
    row2 = [1.2879   -1.6084    0.3205         0         0]
    row3 = [ 0    1.2821   -4.6154    3.3333         0]
    row4 = [0         0   15.0000  -21.6667    6.6667]
    row5 = [0         0         0    8.8889   -8.8889]
    @test μ * L_1_minus + σ^2/2 * L_2 ≈ cat(row1, row2, row3, row4, row5, dims = 1) atol = 1e-4
    row1 = [-1.3223    1.3223         0         0         0]
    row2 = [0.3788   -1.4685    1.0897         0         0]
    row3 = [0    0.5128   -8.8462    8.3333         0]
    row4 = [0         0   10.0000  -20.0000   10.0000]
    row5 = [ 0         0         0    5.5556   -5.5556]
    @test -μ * L_1_plus + σ^2/2 * L_2 ≈ cat(row1, row2, row3, row4, row5, dims = 1) atol = 1e-3 # Test for positive drift
end 