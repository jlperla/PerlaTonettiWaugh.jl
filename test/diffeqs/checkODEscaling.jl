#=
    # Parameters
        const ξ = SOMETHING # Constant rescaling value. 
        const r = SOMETHING # Constant state multiplier
        const μ = SOMETHING # Constant drift < 0. 
        const σ = SOMETHING # Constant shock
        const π(ξ) = z -> exp(ξ * z) # Rescaled π function. This defines a family of π functions, parametrized by ξ. 

    # Grid 
        x = SOMETHING # Some partition of [0, 5].

    # Vanilla (old) solution 
        x, L_1, L_1_plus, L_2 = diffusionoperators(x) # Regular code. 
        L_T = r*I - μ*L_1 - σ^2/2 * L_2 # Construct the aggregate operator. 
        v = L_T \ π(1.0).(x) # Solution to the differential equation. 

    # New solution. 
        x, L_1, L_1_plus, L_2 = rescaled_diffusionoperators(x, ξ)
        L_T = r*I - μ*L_1 - σ^2/2 * L_2 # Construct the aggregate operator. 
        v = L_T \ π(ξ).(x) # Solution to the rescaled differential equation. 

    # Tests 
        # Do the conversion
        # Check for absolute similarity (norms, etc.)
        # Check for greater stability with the rescaled problem. 
=#