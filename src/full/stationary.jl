# Implements the algebraic stationary solution for the full model. Returns the equilibrium quantities (g, Ω, π) determined by equations H.15-H.17.


# Gives us the (full algebraic) stationary solution for a set of params and an initial x. 
function stationary_algebraic(params, init_x = defaultiv(params); kwargs...)
    @assert params.υ > 0 && params.κ > 0 # Parameter validation
    sol = nlsolve(vals -> stationary_algebraic_aux(vals, params), init_x; inplace = false, xtol = -Inf, kwargs...) # distance in x-space isn't reliable because of ContinuousTransformations
    if ~converged(sol)
        throw(sol) # If it fails, throw the results as an exception. 
    end
    g, z_hat, Ω  = sol.zero
    @assert z_hat > 1 && Ω > 0 && g > 0 # Validate solutions. 
    staticvalues = staticvals(sol.zero, params)
    return merge(staticvalues, (g = g, z_hat = z_hat, Ω = Ω))
end 

# Gives us the residuals for a point x in state-space and a set of params. 
function stationary_algebraic_aux(vals, params)    
    # Grab values and intermediate quantities. 
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    g, z_hat, Ω = vals
    @unpack π, π_min, z_bar, L_tilde, S = staticvals(vals, params)
    # Validate parameters. 
    # Compute intermediate algebraic quantities. 
    r = ρ + γ*g + δ # H.6 (WORKING PAPER)
    a = inv(r - g - (σ-1)*(μ-g+(σ-1)*υ^2/2)) # H.4 (WORKING PAPER)
    ν = (μ-g)/υ^2 + sqrt((g-μ)^2/υ^4 + (r-g)/(0.5 * υ^2)) # H.3 (WORKING PAPER)
    b = (1 - a*(r-g))*d^(1-σ)*z_hat^(ν+σ-1) # H.5 (WORKING PAPER)
    # Calculate and assign residuals. 
        big_denom = ν*(θ + ν)*(θ - σ + 1) # Part of H.16 (WORKING PAPER)
        denom_1 = a*(g - r) # Part of H.16 (WORKING PAPER)
        num_1 = ν*(N-1)*(θ - σ + 1)*(d^(1 - σ)*(θ + ν)*z_hat^(-θ + σ - 1)-b*θ*z_hat^(-θ-ν)) # Part of H.16 (WORKING PAPER)
        num_2 = θ*(ν*(N-1)*d^(1-σ)*(θ+ν)*z_hat^(-θ + σ - 1) + (ν + σ - 1)*(θ + ν - σ + 1)) # Part of H.16 (WORKING PAPER)
    return [x/π_min - a*(χ/(1-χ))*(σ + ν - 1)/ν, # H.15 (WORKING PAPER)
    1 + (σ-1)/ν - (num_1/denom_1 + num_2)/big_denom + (χ/(1-χ))*(σ + ν - 1)/(ν), # H.16 (WORKING PAPER)
    π_min - (1- L_tilde)/((σ -1)*z_bar^(σ-1))] # H.17 (WORKING PAPER)
end

function staticvals(vals, params)
    # Unpack x. 
    g, z_hat, Ω = vals
    # Unpack params. 
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    # Compute values. 
    S = θ * (g - μ - θ * υ^2/2) # 28 (PDF)
    L_tilde = Ω*((N-1)*z_hat^(-θ)*κ + ζ*S*δ) # 29 (PDF). Uses E = δ in steady state, from 17
    z_bar = (Ω * (θ/(1 + θ - σ))*(1 + (N-1)*d^(1-σ)*z_hat^(σ-1-θ)))^(1/(σ-1)) # 30 (PDF)
    π_min = (1-L_tilde)/((σ-1)*z_bar^(σ-1)) # 31 (PDF)
    i = z -> (z >= log(z_hat)) # Indicator function used in 32
    π = z -> π_min * (1 + (N-1) * d^(1-σ)*i(z)) - (N-1) * κ * exp(-z * (σ-1)) * i(z) # 32 (PDF)
    return (π = π, π_min = π_min, z_bar = z_bar, L_tilde = L_tilde, S = S)
end 

# Default initial values 
defaultiv(params) = [0.01, 2.0, 1.0]

# Numerical method.
function stationary_numerical(params, z, init_x = defaultiv(params); kwargs...)
    # Unpack params and settings. 
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
    # Validate parameters 
    @assert params.υ > 0 
    @assert params.κ > 0 

    # Discretization objects and quadrature weights. 
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # Operators. 
    ω = ω_weights(z, θ, σ-1) # Get quadrature weights for the distribution on the rescaled grid. 

    # Set up the transformations. 
    map_circle = x -> x/sqrt(1+x*x)
    map_g = x -> (map_circle(x) + 1.1)*(ρ) # Take real values into g values in [1e-10, 0.99ρ]
    map_z_hat = x -> (map_circle(x) + 1.1)*(10) # Fix upper bound at 10 arbitrarily. Since z_bar = f(z_hat), which we're solving for. 
    map_Ω = x -> (map_circle(x) + 1.1)*10 # Takes real values into Ω values in [0, 10]

    # Define the system of equations we're solving.
    function stationary_numerical_given_vals(vals)
        # Unpack arguments and apply ContinuousTransformations
        g_raw, z_hat_raw, Ω_raw = vals 
        g = map_g(g_raw)
        z_hat = map_z_hat(z_hat_raw)
        Ω = map_Ω(Ω_raw) 
        # Compute static numerical values. 
        @unpack π, π_min, z_bar, L_tilde, S = staticvals(vals, params) # Grab static values. 
        # Set the backward-looking derivative. 
        L_tilde_derivative = 0.0 # In the steady state. 
        # Compute more interim quantities. 
        ρ_tilde = ρ + δ + (σ-1)*(μ - g + (σ-1)*υ^2/2) # 33 (PDF)
        A = ρ_tilde*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2 # 23 (PDF)
        # We know from (24, PDF) that A*v = π.(z), because in steady state v'(t) = 0
        v = π.(z) \ A # This is rescaled. 

        #= 
        Return the residuals from the equations to be solved.
        =#

        # Value-matching condition. 
        value_matching = v[1] - dot(v, ω) + ζ # 25 (PDF)
       
        # Export threshold condition. 
        export_threshold = z_hat^(σ-1) - κ*d^(σ-1)*π_min^(-1) # 26 (PDF)

        # Free-entry condition. 
        free_entry = v[1] - ζ*(1-χ)/(χ) # 27 (PDF)
        
        return [value_matching, free_entry, export_threshold]
    end 

    # Solve the equation for the steady state.
    sol = nlsolve(stationary_numerical_given_vals, init_x; inplace = false, xtol = -Inf, kwargs...) # xtol = -Inf because distances in x-space are unreliable, due to ContinuousTransformations. 

    # Error if not converged
    if ~converged(sol)
        throw(sol) # If it fails, throw the results as an exception. 
    end

    # Grab and map the values if converged.
    g_T_raw, z_hat_T_raw, Ω_T_raw = sol.zero
    g_T = map_g(g_T_raw)
    z_hat_T = map_z_hat(z_hat_T_raw)
    Ω_T = map_Ω(Ω_T_raw)

    # Recreate the stationary_numerical_given_vals from those values. 


    return merge(staticvalues, (g = g_T, z_hat = z_hat_T, Ω = Ω_T, v_tilde = v_tilde))
end 