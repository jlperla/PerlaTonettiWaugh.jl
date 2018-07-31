#= 
Test on reascaling method and make comparison
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

    # Vanilla (old) solution 
        x, L_1, L_1_plus, L_2 = diffusionoperators(x) # Regular code. 
        L_T = r*I - μ_val2*L_1 - σ_val2^2/2 * L_2 # Construct the aggregate operator. 
        v = L_T \ π(ξ).(x) # Solution to the differential equation. 

    # New solution. 
        x, L_1_tild, L_1_plus_tild, L_2_tild = rescaled_diffusionoperators(x, ξ)
        # modify r and μ_val2 ,π to match (13)
        r_tild = r - ξ*μ_val2 - σ_val2^2/2*ξ^2;
        μ_tild = μ_val2 + σ_val2^2*ξ;
        π_tild(x) = 1;
        L_T_tild = r_tild*I - μ_tild*L_1_tild - σ_val2^2/2 * L_2_tild # Construct the aggregate operator. 
        v_tild = L_T_tild \ π_tild.(x) # Solution to the rescaled differential equation.        

    # Tests 
        # Do the conversion
        v_rescale = v_tild.* exp.(ξ*x)
        # Check for absolute similarity (norms, etc.)
        @show norm(v-v_rescale)
        @show norm(v-v_rescale,Inf)

        # Check for greater stability with the rescaled problem. 
        x_alt = linspace(x_min,x_max,M+100)
        x_alt, L_1_tild, L_1_plus_tild, L_2_tild = rescaled_diffusionoperators(x_alt, ξ)
        # modify r and μ_val2 ,π to match (13)
        r_tild = r - ξ*μ_val2 - σ_val2^2/2*ξ^2;
        μ_tild = μ_val2 + σ_val2^2*ξ;
        π_tild(x) = 1;
        L_T_tild = r_tild*I - μ_tild*L_1_tild - σ_val2^2/2 * L_2_tild # Construct the aggregate operator. 
        v_tild_alt = L_T_tild \ π_tild.(x_alt) # Solution to the rescaled differential equation. 
        
        # test the pre-conversion difference of v_tild
        v_int_tild=LinInterp(x, v_tild)
        @show norm(v_tild_alt-v_int_tild.(x_alt))

        # conversion to v
        v_rescale_alt = v_tild_alt.*exp.(ξ*x_alt)
        v_int=LinInterp(x, v_rescale)
        @show norm(v_rescale_alt-v_int.(x_alt))

        # use the ordinary method
        x_alt, L_1, L_1_plus, L_2 = diffusionoperators(x_alt) # Regular code. 
        L_T = r*I - μ_val2*L_1 - σ_val2^2/2 * L_2 # Construct the aggregate operator. 
        v_alt = L_T \ π(ξ).(x_alt) # Solution to the differential equation. 
        v_int_old=LinInterp(x, v)
        @show norm(v_alt-v_int_old.(x_alt))

        # Plot of v 
        plot(x,v_tild)
        plot!(x_alt,v_tild_alt)