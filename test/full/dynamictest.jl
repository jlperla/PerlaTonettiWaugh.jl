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

# regression tests, based on version 6453342
@testset "Regression tests for f! in `solve_dynamics`" begin
    rng_seed = MersenneTwister(123456) # set RNG seed for reproducibility
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
        @test resid[1] ≈ -0.036653936375744864
        @test resid[2] ≈ 0.02145051714277048
        @test resid[10] ≈ 0.0004069292969686722
        @test resid[M] ≈ -4.439559001072346e-5
        @test resid[M+1] ≈ -3.1764554464519534e-5
        @test resid[M+2] ≈ -0.010535368447424975
    end

    @testset "t = T / 2" begin
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
        @test resid[1] ≈  -0.04387774636852337
        @test resid[2] ≈  0.022224588053667838
        @test resid[10] ≈ -0.00016380102265192709
        @test resid[M] ≈ -0.00022529945836316512
        @test resid[M+1] ≈ -2.802589158656943e-5
        @test resid[M+2] ≈  0.008324583711737166
    end
    
    @testset "t = 0.1" begin
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
        @test resid[1] ≈ 0.09009608302866695
        @test resid[2] ≈ 0.24485202918954876
        @test resid[10] ≈ -0.0006604034953704707
        @test resid[M] ≈ -0.0005796653180272643
        @test resid[M+1] ≈ 0.006699045899212219
        @test resid[M+2] ≈ -0.001569813559400668
    end


end

# Solve and compute residuals
@time solved = solve_dynamics(params_T, stationary_sol_T, settings, T, Ω)

@test mean(mean(solved.residuals[:,1:M], dims = 1)) ≈ 0 atol = 1e-03 # mean residuals for system of ODEs
@test mean(mean(solved.residuals[:,(M+1)])) ≈ 0 atol = 1e-03 # mean residuals for value matching condition
@test mean(mean(solved.residuals[:,(M+2)])) ≈ 0 atol = 1e-03 # mean residuals for export threshold condition

v_hat_t0 = map(z -> exp((σ-1)*z), z_grid) .* solved.v[1]
@test any((v -> v < 0).(diff(v_hat_t0))) == false # after reparametrization v_hat should be increasing