# TODO: ADD `solve_ptw_full` to perform optimization for Ω
function entry_residuals(params_T, stationary_sol_T, settings, T, Ω_vec0, Ω_nodes, entry_residuals_nodes)
    @unpack δ, χ, ζ = params_T
    @unpack z = settings
    M = length(z)
    @assert length(entry_residuals_nodes) >= length(Ω_nodes) - 2 # there should be enough sample points
    @assert Ω_nodes[1] ≈ 0
    @assert Ω_nodes[end] ≈ T


    entry_residuals = similar(entry_residuals_nodes) # entry_residuals from equil. Ω, evaluated at entry_residuals_nodes

    ## TODO: add optimization step; each candidate for solution is denoted as Ω_vec.  
    ## find the equilibrium Ω based on entry_residuals
    Ω_vec = Ω_vec0 # TODO: remove this -- this assumes that the first guess for Ω_vec is used throughout the code
    solved_dynamics = solve_dynamics_by_vector_Ω(params_T, stationary_sol_T, settings, T, Ω_vec, Ω_nodes)

    t = solved_dynamics.sol.t # FIXME: later make solved_dynamics return t so that we don't need to call the entire sol
    v = solved_dynamics.v
    if (!issorted(t)) # FIXME: later make solved_dynamics return t in sorted way
        # this assumes that t and v are in the same orders 
        v = sort(sortperm(t))
        t = sort(t)
    end
    v_interpolation = LinInterp(t, v)
    # compute entry residuals from the solution
    entry_residuals_vec = map(t -> v_interpolation(t)[1] - ζ * (1-χ) / χ, entry_residuals_nodes)
    ## end of optimization

    # perform linear interpolation on entry_residuals 
    entry_residuals_interpolation = LinInterp(entry_residuals_nodes, entry_residuals_vec)

    return (entry_residuals_interpolation = entry_residuals_interpolation,
            entry_residuals = entry_residuals_vec, solved_dynamics = solved_dynamics)
end

function solve_dynamics_by_vector_Ω(params_T, stationary_sol_T, settings, T, Ω_vec, Ω_nodes)
    # interpolate Ω based on candidate Ω_vec
    Ω_interpolation_instance = LinInterp(Ω_nodes, Ω_vec) # perform linear interpolation
    Ω_interpolation(t) = Ω_interpolation_instance(t) # return interpolated Ω based on Ω_vec

    # solve ptw problem based on the interpolated Ω
    return solve_dynamics(params_T, stationary_sol_T, settings, T, Ω_interpolation)
end

function solve_dynamics(params_T, stationary_sol_T, settings, T, Ω)
    @unpack δ, N, σ, θ, d = params_T
    @unpack z, tstops, Δ_E = settings
    @assert params_T.γ ≈ 1 # γ has to be close 1 to have consistent results with the stationary solutions
    M = length(z)

    # define E(t) based on FD
    E(t) = (log(Ω(t+Δ_E)) - log(Ω(t-Δ_E)))/(2*Δ_E) + δ

    # define the corresponding DAE problem
    p = get_p(params_T, stationary_sol_T, settings, Ω, T)
    dae = PTW_DAEProblem(params_T, stationary_sol_T, settings, E, Ω, T, p)

    # solve solutions
    sol = DifferentialEquations.solve(dae.dae_prob, callback = dae.cb) # solve! # TODO: include tstops = tstops in the argument.
    @unpack u, du, t = sol

    residuals = zeros(length(t), length(u[1]))
    equilibriums = []
    p_for_extraction_template = (; [pair for pair in pairs(p) if pair[1] != :ts && pair[1] != :vals]...) # copy, removing ts and vals 
    for (i, t) in enumerate(t)
        # change p so that the last element in ts/vals is from one time step forward, not the end of the entire history 
        p_for_extraction = merge(p_for_extraction_template, 
                                (ts = [p.ts[max(1, (i-1))]],
                                vals = [p.vals[max(1, i-1)]]))

        # compute residual at t
        residual = zeros(length(u[1]))
        dae.dae_prob.f(residual, du[i], u[i], p_for_extraction, t)
        residuals[i,:] = residual
        # compute stationary equilibrium at t
        v_1 = u[i][1]
        g = u[i][M+1]
        z_hat = u[i][M+2]
        # compute λ_ii and c
        equilibrium = dae.stationary_equilibrium(v_1, g, z_hat, E(t), Ω(t), t)
        λ_ii = 1 / (1 + (N-1)*z_hat^(σ-1-θ)*d^(1-σ))
        c = (θ / (1-σ+θ))^(1/(σ-1))*(1-equilibrium.L_tilde)*Ω(t)^(1/(σ-1))*λ_ii^(1/(1-σ))
        # TODO: add log_M and U later
        equilibrium = merge(equilibrium, (λ_ii = λ_ii, c = c, E = E(t),))
        push!(equilibriums, equilibrium)
    end

    # extract solutions
    v = map(u -> u[1:M], u)
    g = map(u -> u[M+1], u)
    z_hat = map(u -> u[M+2], u)
    S = map(eq -> eq.S, equilibriums)
    L_tilde = map(eq -> eq.L_tilde, equilibriums)
    z_bar = map(eq -> eq.z_bar, equilibriums)
    π_min = map(eq -> eq.π_min, equilibriums)
    π_tilde = map(eq -> eq.π_tilde, equilibriums)
    λ_ii = map(eq -> eq.λ_ii, equilibriums)
    c = map(eq -> eq.c, equilibriums)
    E = map(eq -> eq.E, equilibriums)
    entry_residual = map(eq -> eq.entry_residual, equilibriums)

    return (v = v, g = g, z_hat = z_hat,
            S = S, L_tilde = L_tilde, z_bar = z_bar, π_min = π_min, π_tilde = π_tilde,
            λ_ii = λ_ii, c = c, E = E, entry_residual = entry_residual,
            t = t, p = p, sol = sol, f! = dae.dae_prob.f, residuals = residuals, equilibriums = equilibriums)
end

# Implementation of the full model with time-varying objects, represented by DAE
function PTW_DAEProblem(params_T, stationary_sol_T, settings, E, Ω, T, p)
    # Unpack params and settings.
    @unpack z = settings
    M = length(z)
    @unpack L_1, L_2, z, M, T, μ, υ, σ, d, κ, ω, θ, δ, χ, N, ζ, ρ = p

    get_S(g) = θ * (g - μ - θ * υ^2/2)
    get_L_tilde(S, z_hat, E, Ω) = Ω * ((N-1) * z_hat^(-θ)*κ + ζ*(S + E / χ))

    function stationary_equilibrium(v_1, g, z_hat, E, Ω, t)
        S = get_S(g)
        L_tilde = get_L_tilde(S, z_hat, E, Ω)
        z_bar = Ω * (θ / (1 + θ - σ)) * (1 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ))
        π_min = (1 - L_tilde) / ((σ-1)*z_bar)
        π_tilde = π_min * (1.0.+(N-1)*d^(1-σ)*(z .>= log(z_hat))) - (N-1)*κ*exp.(-(σ-1).*z).*(z .>= log(z_hat))
        entry_residual = v_1 - ζ * (1-χ) / χ
        return (S = S, L_tilde = L_tilde, z_bar = z_bar,
                π_min = π_min, π_tilde = π_tilde,
                entry_residual = entry_residual)
    end

    cb = FunctionCallingCallback((u, t, integrator) -> (push!(p.ts, t); push!(p.vals, get_L_tilde(get_S(u[M+1]), u[M+2], E(t), Ω(t)))), tdir = -1, func_start = false)

    # Dynamic calculations, defined for each time ∈ t.
    function f!(residual,du,u,p,t)
        residual .= 0

        # Carry out calculations.
        # v = u[1:M]; note that this line is disabled and `v` is replaced with u[1:M] to prevent huge memory allocation 
        g = u[M+1]
        z_hat = u[M+2]
        x = ζ
        Ω_t = Ω(t)
        E_t = E(t)

        S = get_S(g)
        L_tilde = get_L_tilde(S, z_hat, E_t, Ω_t)
        z_bar = Ω(t) * (θ / (1.0 + θ - σ)) * (1.0 + (N-1) * d^(1-σ) * z_hat^(σ-1-θ))
        π_min = (1 - L_tilde) / ((σ-1)*z_bar)
        π_tilde = π_min * (1.0.+(N-1)*d^(1-σ)*(z .>= log(z_hat))) - (N-1)*κ*exp.(-(σ-1).*z).*(z .>= log(z_hat))

        # compute the derivative of L_tilde
        L_tilde_derivative_term = 0
        if (t < T - eps())
            t_forward = p.ts[end]
            L_tilde_forward = p.vals[end]
            L_tilde_derivative_term = (log(1 - L_tilde_forward) - log(1-L_tilde))/(t_forward-t) # Reverse direction. 
        end

        # Form the DAE at t.
        # Note that
        # A_t = (ρ + δ + L_tilde_derivative_term - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2
        # and we are decomposing this into the terms involving (I with L_1) AND (L_2) -- otherwise it will perform elementwise operation 
        # on two banded matrices which takes extreme amount of memory that are not needed.
        residual[1:M] = (ρ + δ + L_tilde_derivative_term - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*u[1:M] # system of ODEs (eq:28)
        residual[1:M] .-= (μ - g + (σ-1)*υ^2)*L_1*u[1:M] 
        residual[1:M] .-= (υ^2/2)*L_2*u[1:M] 
        residual[1:M] .-= du[1:M]
        residual[1:M] .-= π_tilde
        residual[M+1] = u[1] + x - dot(ω, u[1:M]) # residual (eq:25)
        residual[M+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min # export threshold (eq:31)
    end

    u0 = [p.v_T; p.g_T; p.z_hat_T]
    du0 = zeros(M+2)

    return (dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), differential_vars = [trues(M); false; false], p),
            cb = cb, stationary_equilibrium = stationary_equilibrium, f! = f!)
end

# return the parameters and functions needed to define dynamics
function get_p(params_T, stationary_sol_T, settings, Ω, T)
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = params_T
    @unpack z = settings
    M = length(z)

    # Unpack the stationary solution.
    v_T = stationary_sol_T.v_tilde
    g_T = stationary_sol_T.g
    z_hat_T = stationary_sol_T.z_hat
    Ω_T = stationary_sol_T.Ω

    # Quadrature weighting
    ω = ω_weights(z, θ, σ-1)  # TODO: compute integral for eq 30

    # Discretize the operator.
    z, L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z, σ-1) # L_1_minus ≡ L_1 is the only one we use.

    # Pre-load the arrays with the stationary values. 
    L_tilde_T = stationary_sol_T.L_tilde 

    # Create the arrays to hold ts and L_tildes. 
    ts = Array{Float64}([T])
    vals = Array{Float64}([L_tilde_T])
        
    # Bundle as before.
    p = (L_1 = L_1_minus, L_2 = L_2, z = z, N = N, M = M, T = T, θ = θ, σ = σ, κ = κ,
        ζ = ζ, d = d, ρ = ρ, δ = δ, μ = μ, υ = υ, χ = χ, ω = ω, Ω = Ω,
        v_T = v_T, g_T = g_T, z_hat_T = z_hat_T, Ω_T = Ω_T,
        ts = ts, vals = vals) #Named tuple for parameters.
    return p
end
