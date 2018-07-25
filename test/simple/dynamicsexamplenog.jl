# Global constants. 
    solve(x,y)=Sundials.solve(x,y)
    # For the test grid. 
    x_min = 0.01
    x_max = 1.0
    M = 20
    x=linspace(x_min,x_max,M) # Since we only use the grid now. 
    # For the time boundary. 
    T = 10.0
    # Individual model parameters. 
    rho = 0.15
    sigma_bar = 0.1
    c_tilde(t, x) = x + 0.0001*t
    sigma_tilde(t, x) =  sigma_bar
    mu_tilde(t, x) = 0.1*x *(1.0 + 4.0 * t)
    # Solver settings. 
    basealgorithm = CVODE_BDF() 
    # plotevery = 5 irrelevant now, since there is no more plotting. 

#= 
    Tests on vanilla examples. 
=#
    # Vanilla example ODE. 
    prob = createsimpleODEproblem(c_tilde, sigma_tilde, mu_tilde, x,  T, rho)
    sol = solve(prob, basealgorithm)
    @test issorted(sol[end])

    # Vanilla example DAE.
    probDAE = createsimpleDAEproblem(c_tilde, sigma_tilde, mu_tilde, x,  T, rho)
    solDAE = solve(probDAE, IDA())
    @show(issorted(solDAE[end][1:M]))

    #Check they are "reasonably" close
    @show norm(sol[1] - solDAE[1][1:M])
    @show norm(sol[end] - solDAE[end][1:M])

#=
    Test for nonuniform. 
=#
    
    # Create new grid. 
    x_add=linspace(x_min,x_max,10)
    x_comb=unique(sort([x;x_add]))
    M_comb = length(x_comb)

    # Test new ODE. 
    prob_nonuni = createsimpleODEproblem(c_tilde, sigma_tilde, mu_tilde, x_comb, T, rho)
    sol_nonuni = solve(prob_nonuni, basealgorithm)
    @test issorted(sol_nonuni[end])

    # Test its closeness to the uniform ODE. 
    @show norm(sol[1,1]-sol_nonuni[1,1])
    @show norm(sol[end,end]-sol_nonuni[end,end])

    # interpolate uniform onto non uniform grid 
    sol_int=LinInterp(x, sol[end])
    @show norm(sol_int.(x_comb)-sol_nonuni[end])

    # Test new DAE. 
    probDAE_nonuni = createsimpleDAEproblem(c_tilde, sigma_tilde, mu_tilde, x_comb, T, rho)
    solDAE_nonuni = solve(probDAE_nonuni, IDA())
    @show(issorted(solDAE_nonuni[end][1:M_comb]))

    #Check they are "reasonably" close
    @show norm(sol_nonuni[1] - solDAE_nonuni[1][1:M_comb])
    @show norm(sol_nonuni[end] - solDAE_nonuni[end][1:M_comb])

#= 
    Backwards tests. 
=#
    # Test for ODE with backwards drift. 
    mu_tilde(t, x) = -1 * 0.1*x *(1.0 + 4.0 * t)
    probBackODE = createsimpleODEproblem(c_tilde, sigma_tilde, mu_tilde, x,  T, rho)
    solBackODE = solve(probBackODE, basealgorithm)

    # Monotonicity test
    @test issorted(solBackODE[end])
    
    # Some invariance/regression tests on the backwards solution
    @test solBackODE.u[1] ≈ [0.17211, 0.180156, 0.189957, 0.200501, 0.211449, 0.222651, 0.234025, 0.245525, 0.257119, 0.268787, 0.280515, 0.292291, 0.304109, 0.315961, 0.327842, 0.339749, 0.351679, 0.363628, 0.375588, 0.387297] atol = 1e-4
    @test solBackODE.u[5] ≈ [0.172132, 0.180285, 0.190226, 0.200922, 0.212028, 0.223389, 0.234924, 0.246585, 0.258342, 0.270173, 0.282063, 0.294003, 0.305983, 0.317998, 0.330043, 0.342113, 0.354206, 0.366318, 0.378442, 0.390302] atol = 1e-4
    @test solBackODE.alg == Sundials.CVODE_BDF{:Newton,:Dense}(0, 0, 0, 0, false, 10,
    5, 7, 3, 10)