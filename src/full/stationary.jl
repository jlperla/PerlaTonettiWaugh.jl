# Default initial values 
defaultiv(params) = [0.01, 2.0, 1.0]

# Gives us the (full algebraic) stationary solution for a set of params and an initial x. 
function stationary_algebraic(params, init_x = defaultiv(params); kwargs...)
  # Unpack parameters. 
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params

  # Validate parameters. 
    @assert υ > 0 
    @assert κ > 0 
    @assert γ ≈ 1 # For the same reason we do this in the dynamics. 

  # Define the static equilibrium objects. It's easiest to do this as a closure, because of the dependencies between the equations.
    function staticvals(vals)
      # Unpack argument. 
      g, z_hat, Ω = vals
      # Compute interim quantities.
      F(z) = 1 - z^(-θ) # H.1 
      r = ρ + γ*g + δ # H.6 
      ν = (μ-g)/υ^2 + sqrt(((g-μ)/υ^2)^2 + (r-g)/(υ^2/2)) # H.3
      a = (r - g - (σ - 1)*(μ - g + (σ - 1)*υ^2/2))^(-1) # H.4 
      b = (1 - a*(r-g))*d^(1-σ)*z_hat^(ν + σ - 1) # H.5 
      S = θ * (g - μ - θ * υ^2 /2) # H.2   
      L_tilde = Ω * ((N-1)*(1-F(z_hat))*κ + (1-η)*ζ*(S + δ/χ)) # H.7
      z_bar = (Ω * (θ/(1 + θ - σ) + (N-1)*(1-F(z_hat))*d^(1-σ)*(z_hat^(-1 + σ)*θ/(1 + θ - σ))))^((σ-1)^(-1)) # H.8
      w = σ^(-1)*z_bar # H.10
      x = ζ * (1- η + η * Theta / w) # H.11
      π_min = (d^(σ-1) * κ)/(z_hat^(σ-1)) # Inversion of H.9    
      return (F = F, r = r, ν = ν, a = a, b = b, S = S, L_tilde = L_tilde, z_bar = z_bar, w = w, x = x, π_min = π_min)
    end

  # Define the function we'll pass to the solver. 
    function stationary_algebraic_given_vals(vals)
      # Grab the static values. 
        @unpack F, r, ν, a, b, S, L_tilde, z_bar, w, x, π_min = staticvals(vals)
      # Unpack the argument. 
        g, z_hat, Ω = vals
      # Validate arguments. 
        @assert z_hat > 1 
        @assert Ω > 0 
        @assert g > 0 
      # Calculate and assign residuals. 
        big_denom = ν*(θ + ν)*(θ - σ + 1) # Part of H.16
        denom_1 = a*(g - r) # Part of H.16
        num_1 = ν*(N-1)*(θ - σ + 1)*(d^(1 - σ)*(θ + ν)*z_hat^(-θ + σ - 1)-b*θ*z_hat^(-θ-ν)) # Part of H.16
        num_2 = θ*(ν*(N-1)*d^(1-σ)*(θ+ν)*z_hat^(-θ + σ -1) + (ν + σ - 1)*(θ + ν - σ + 1)) # Part of H.16
        return [x/π_min - a*(χ/(1-χ))*(σ + ν - 1)/ν, 1 + (σ-1)/ν - (num_1/denom_1 + num_2)/big_denom + (χ/(1-χ))*(σ + ν - 1)/(ν), π_min - (1- L_tilde)/((σ -1)*z_bar^(σ-1))]
    end 

    # Run the solver. 
      sol = nlsolve(vals -> stationary_algebraic_given_vals(vals, params), init_x; inplace = false, kwargs...)

    # Throw error if we didn't converge. 
    if ~converged(sol)
        throw(sol) # If it fails, throw the results as an exception. 
    end
    
    # Return results if we did. 
      g, z_hat, Ω  = sol.zero
      @assert z_hat > 1 
      @assert Ω > 0 
      @assert g > 0 # Final results check. 
      staticvalues = staticvals(sol.zero, params)
      return merge(staticvalues, (g = g, z_hat = z_hat, Ω = Ω))
end 

# Numerical method.
function stationary_numerical(params, z, init_x = defaultiv(params); kwargs...)
    # Unpack arguments. 
      @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params
      @assert params.υ > 0 && params.κ > 0 # Parameter validation

    # Validate arguments. 
      @assert υ > 0 
      @assert κ > 0 
      @assert γ ≈ 1 # For the same reason we do this in the dynamics. 

    # Discretization objects and quadrature weights. 
      z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # Operators. 
      L_1 = L_1_minus # L_1 ≡ L_1_minus. 
      ω = ω_weights(z, θ, σ-1) # Get quadrature weights for the distribution on the rescaled grid. 

    # Set up the continuous transformations. 
      map_circle = x -> x/sqrt(1+x*x)
      map_g = x -> (map_circle(x) + 1 + 1e-10)*(ρ) # Take real values into g values in [1e-10, 0.99ρ]
      map_z_hat = x -> (map_circle(x) + 1 + 1e-5)*(10) # Fix upper bound at 10 arbitrarily. Since z_bar = f(z_hat), which we're solving for. 
      map_Ω = x -> (map_circle(x) + 1)*10 # Takes real values into Ω values in [0, 10]

    # Set up static value calculations. 
      function staticvals(vals)
          # Unpack arguments. 
            g, z_hat, Ω = vals 
          # Start calculations. 
            S = θ*(g - μ - θ*(υ^2/2)) # 28 in the PDF 
            E = δ # 17 in the PDF, at steady state 
            L_tilde = Ω*((N-1)*z_hat^(-θ)*κ + ζ*(S + E/χ)) # 29 in the PDF 
            z_bar = (Ω*(θ/(1 + θ - σ))*(1 + (N-1)*d^(1-σ)*z_hat^(σ-1-θ)))^(1/(σ-1)) # 30 in the PDF 
            π_min = (1-L_tilde)/((σ-1)*z_bar^(σ-1)) # 31 in the PDF 
            i = z -> (z >= log(z_hat)) # Indicator function from 32 in the PDF 
            π = z -> π_min*(1 + (N-1)*d^(1-σ))*i(z) - (N-1)*κ*exp(-(σ-1)*z)*i(z) # 32 in the PDF 
            ρ_tilde = ρ + δ - (σ-1)*(μ-g+(σ-1)*(υ^2/2)) # 33 in the PDF, at steady state. 
          # Return. 
            return (S = S, E = E, L_tilde = L_tilde, z_bar = z_bar, π_min = π_min, π = π, ρ_tilde = ρ_tilde)
      end 

    # Define the system of equations we're solving.
      function stationary_numerical_given_vals(vals)
        # Unpack and transform the raw arguments. 
          g_raw, z_hat_raw, Ω_raw = vals 
          g = map_g(g_raw)
          z_hat = map_z_hat(z_hat_raw)
          Ω = map_Ω(Ω_raw)
        # Validate the variables. 
          @assert g > 0 
          @assert z_hat > 1
          @assert Ω > 0 
        # Grab the static values. 
          @unpack S, E, L_tilde, z_bar, π_min, π, ρ_tilde = staticvals([g, z_hat, Ω])
        # Define the discretized aggregate operator. 
          A = ρ_tilde * I - (μ - g + (σ-1)*υ^2)*L_1 - υ^2/2 * L_2 # 23 in the PDF
        # Calculate the value function. 
          v = A \ π.(z) # 24 in the PDF, at steady state
        # Apply the conditions.
          value_matching = v[1] - dot(v, ω) + ζ # 25 in the PDF (VALUE MATCHING CONDITION)
          export_threshold = z_hat^(σ-1) - κ*d^(σ-1)*π_min^(-1) # 26 in the PDF (EXPORT THRESHOLD CONDITION)
          free_entry = v[1] - ζ*(1-χ)/χ # 27 in the PDF, at steady state (E = δ > 0)
        # Return. 
          return [value_matching, export_threshold, free_entry]
      end

    # Solve the above equation. 
      sol = nlsolve(stationary_numerical_given_vals, init_x; inplace = false, xtol = -Inf, kwargs...)
    
    # Throw an error if we failed to converge. 
      if ~converged(sol)
        throw(sol) # Will throw all the data in sol as an error. 
      end 

    # Recreate the economic calculations from the solver results. 
      g_T_raw, z_hat_T_raw, Ω_T_raw = sol.zero # Grab raw values. 
      g_T = map(g_T_raw) # Transform raw values. 
      z_hat_T = map(z_hat_T_raw)
      Ω_T = map(Ω_T_raw)
    # Validate parameters as a final check. 
      @assert g_T > 0 
      @assert z_hat_T > 1 
      @assert Ω_T > 0 
    # Grab the static values.
      @unpack S, E, L_tilde, z_bar, π_min, π, ρ_tilde = staticvals([g_t, z_hat_T, Ω_T])
    # Define the discretized aggregate operator. 
      A_T = ρ_tilde * I - (μ - g_T + (σ-1)*υ^2)*L_1 - υ^2/2 * L_2 # 23 in the PDF
    # Calculate the value function. 
      v_T = A_T \ π.(z) # 24 in the PDF, at steady state
    # Apply the conditions.
      value_matching = v_T[1] - dot(v_T, ω) + ζ # 25 in the PDF (VALUE MATCHING CONDITION)
      export_threshold = z_hat_T^(σ-1) - κ*d^(σ-1)*π_min^(-1) # 26 in the PDF (EXPORT THRESHOLD CONDITION)
      free_entry = v_T[1] - ζ*(1-χ)/χ # 27 in the PDF, at steady state (E = δ > 0)
    # Return. 
      return (g = g_T, z_hat = z_hat_T, Ω = Ω_T, v = v_T, S = S, L_tilde = L_tilde, π_min = π_min)
end 