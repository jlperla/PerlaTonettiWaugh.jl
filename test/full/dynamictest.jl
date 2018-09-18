# State grid. 
z_min = 0.0 
z_max = 5.0
M = 1000
z_grid = range(z_min, stop = z_max, length = M) # Since we only care about the grid. 

# Define common objects. 
d_0 = 5
d_T = 2.3701
params = (ρ = 0.02, σ = 4.2508, N = 10, θ = 5.1269, γ = 1.00, κ = 0.013, ζ = 1, η = 0, Theta = 1, χ = 1/(2.1868), υ = 0.0593, μ = 0, δ = 0.053) # Baselines per Jesse. 
σ = params.σ
δ = params.δ

# Compute the stationary solution at t = 0 and t = T first
params_0 = merge(params, (d = d_0,)) # parameters to be used at t = 0
params_T = merge(params, (d = d_T,)) # parameters to be used at t = T

stationary_sol_0 = stationary_numerical(params_0, z_grid) # solution at t = 0
Ω_0 = stationary_sol_0.Ω

stationary_sol_T = stationary_numerical(params_T, z_grid) # solution at t = T
v_T = stationary_sol_T.v_tilde
g_T = stationary_sol_T.g
z_hat_T = stationary_sol_T.z_hat
Ω_T = stationary_sol_T.Ω

# compute the resulting end time and function of Ω
T = 20.0
Ω(t) = Ω_T

settings = (z = z_grid, tstops = 0:1e-3:T, Δ_E = 1e-06)
E(t) = (log(Ω(t+settings.Δ_E)) - log(Ω(t-settings.Δ_E)))/(2*settings.Δ_E) + params.δ
p = PerlaTonettiWaugh.get_p(params_T, stationary_sol_T, settings, Ω, T);
dae = PerlaTonettiWaugh.PTW_DAEProblem(params_0, stationary_sol_T, settings, E, Ω, T, p);

# regression tests, based on version b20c067
@testset "Regression tests for f! in `solve_dynamics`" begin
    RNG_SEED_ORIGIN = 123456 # set RNG seed for reproducibility
    @testset "t = T" begin
        u = [stationary_sol_T.v_tilde; stationary_sol_T.g; stationary_sol_T.z_hat]
        du = zeros(M+2)
        resid = similar(u)
        t = T

        # compute residuals
        dae.f!(resid,du,u,p,t)

        # test if residuals are small enough
        @test mean(resid[1:M]) ≈ 0 atol = 1e-8
        @test mean(resid[M+1]) ≈ 0 atol = 1e-8
        @test mean(resid[M+2]) ≈ 0 atol = 1e-8

        # check if values are close enough as before (regression tests)
        @test resid[1] ≈ -9.016777602344206e-11
        @test resid[2] ≈ -9.016033752917707e-11
        @test resid[10] ≈ -9.013649548972325e-11
        @test resid[M] ≈ -1.392167769953545e-10
        @test resid[M+1] ≈ -2.7400304247748863e-13
        @test resid[M+2] ≈ 4.208980275421936e-9

        @test u[1] ≈ 1.186800000000357
        @test u[2] ≈ 1.168072362393182
        @test u[10] ≈ 1.0372263636598893
        @test u[M] ≈ 0.675633431057062
        @test u[M+1] ≈ 0.020419880554626162
        @test u[M+2] ≈ 1.425712535487968
    end

    @testset "t = T - 1e-3" begin
        rng_seed = MersenneTwister(RNG_SEED_ORIGIN)
        u = [stationary_sol_T.v_tilde; stationary_sol_T.g; stationary_sol_T.z_hat]
        du = zeros(M+2)
        resid = similar(u)
        t = T - 1e-3

        # give some changes
        u[2:3] = u[2:3] .+ rand(rng_seed, 2) * 1e-3
        u[M+1:M+2] = u[M+1:M+2] .+ rand(rng_seed, 2) * 1e-3
        du = rand(rng_seed, length(du)) * 1e-3

        # compute residuals
        dae.f!(resid,du,u,p,t)

        # check if values are close enough as before (regression tests)
        @test resid[1] ≈ 4.542052507407285
        @test resid[2] ≈ 4.529892925281381
        @test resid[10] ≈ 4.0020542455673995
        @test resid[M] ≈ 2.6065676236341866
        @test resid[M+1] ≈ -3.1764554464519534e-5
        @test resid[M+2] ≈ -0.010535368447424975
    end

    @testset "t = T / 2" begin
        rng_seed = MersenneTwister(RNG_SEED_ORIGIN)
        u = [stationary_sol_T.v_tilde; stationary_sol_T.g; stationary_sol_T.z_hat]
        du = zeros(M+2)
        resid = similar(u)
        t = T / 2

        # give some changes
        u[2:3] = u[2:3] .+ rand(rng_seed, 2) * 1e-3
        u[M+1:M+2] = u[M+1:M+2] .+ rand(rng_seed, 2) * 1e-3
        du = rand(rng_seed, length(du)) * 1e-3
        du[1] = du[1] + 1e-2

        # compute residuals
        dae.f!(resid,du,u,p,t)

        # check if values are close enough as before (regression tests)
        @test resid[1] ≈ -0.046196065731366004
        @test resid[2] ≈ 0.02190136138358488
        @test resid[10] ≈ 0.0008070940285961892
        @test resid[M] ≈ 0.00021626561191201432
        @test resid[M+1] ≈ -3.1764554464519534e-5
        @test resid[M+2] ≈ -0.010535368447424975
    end
    
    @testset "t = 0.1" begin
        rng_seed = MersenneTwister(RNG_SEED_ORIGIN)
        u = [stationary_sol_T.v_tilde; stationary_sol_T.g; stationary_sol_T.z_hat]
        du = zeros(M+2)
        resid = similar(u)
        t = 0.1

        # give some changes
        u[1:5] = u[1:5] .+ rand(rng_seed, 5) * 1e-2
        u[M+1:M+2] = u[M+1:M+2] .+ rand(rng_seed, 2) * 1e-3
        du = rand(rng_seed, length(du)) * 1e-3
        du[3] = du[3] - 1e-2

        # compute residuals
        dae.f!(resid,du,u,p,t)

        # check if values are close enough as before (regression tests)
        @test resid[1] ≈ -0.16437057505724734
        @test resid[2] ≈ 0.27593757590866175
        @test resid[10] ≈ 0.0005903931927783018
        @test resid[M] ≈  6.730926950911098e-5
        @test resid[M+1] ≈ 0.004653937079329484
        @test resid[M+2] ≈ -0.00815101908664495
    end
end

# regression tests, based on version b20c067
@testset "Regression tests for `solve_dynamics_by_vector_Ω`" begin
    # Solve and compute residuals, now using vectorized Ω
    Ω_nodes = 0:1e-1:T
    Ω_vec = map(t -> Ω(t), Ω_nodes)
    @time solved = PerlaTonettiWaugh.solve_dynamics_by_vector_Ω(params_T, stationary_sol_T, settings, T, Ω_vec, Ω_nodes)

    # solved residuals should be close to zero, but more generous criteria due to discreteness of Ω_vec
    @test mean(mean(solved.residuals[:,1:M], dims = 1)) ≈ 0 atol = 1e-03 # mean residuals for system of ODEs
    @test mean(mean(solved.residuals[:,(M+1)])) ≈ 0 atol = 1e-03 # mean residuals for value matching condition
    @test mean(mean(solved.residuals[:,(M+2)])) ≈ 0 atol = 1e-03 # mean residuals for export threshold condition

    # check if values are close enough as before (regression tests)
    @test solved.residuals[1,1] ≈ -9.016777602344206e-11
    @test solved.residuals[2,2] ≈ 8.515632643479876e-12
    @test solved.residuals[3,M] ≈ -4.96713781217295e-12
    @test solved.residuals[4,M+1] ≈ 1.1857181902996672e-13
end

# regression tests, based on version b20c067
@testset "Regression tests for `entry_residuals`" begin
    # Solve and compute residuals, now using vectorized Ω
    Ω_nodes = 0:1e-1:T
    Ω_vec0 = map(t -> Ω(t), Ω_nodes)
    entry_residuals_nodes = 1e-1:1e-1:(T-1e-1)
    entry_residuals_nodes_last = entry_residuals_nodes[end]
    @time solved = entry_residuals(params_T, stationary_sol_T, settings, T, Ω_vec0, Ω_nodes, entry_residuals_nodes)

    # check if values are close enough as before (regression tests)
    @test solved.entry_residuals[1] ≈ -0.1765540983606555
    @test solved.entry_residuals[2] ≈  -0.16630819672131136
    @test solved.entry_residuals[3] ≈ -0.156062295081967
    @test solved.entry_residuals_interpolation(entry_residuals_nodes[1] + 1e-3) ≈ -0.17645163934426206
    @test solved.entry_residuals_interpolation(entry_residuals_nodes_last) ≈ 7.563199999999956
    @test solved.entry_residuals_interpolation(entry_residuals_nodes_last) ≈ solved.entry_residuals[end]
end