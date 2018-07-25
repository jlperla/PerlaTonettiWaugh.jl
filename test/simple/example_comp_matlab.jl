# Global constants.
    solve(x,y)=Sundials.solve(x,y)
    # For the test grid. 
    x_min = 0.01
    x_max = 1.0
    M = 20
    x = linspaxce(x_min, x_max, M) # Since we use only the grid now. 
    # For the time boundary. 
    T = 10.0
    # Individual model parameters. 
    rho = 0.0286
    sigma_bar = 0.02
    c_tilde(t, x) = x + 0.0*t # All those zeros unnecessary; 0.0 does the job fine. 
    sigma_tilde(t, x) =  sigma_bar
    mu_tilde(t, x) = 0.0*x *(1.0 + 0.0 * t)-0.0161
    # Solver settings. 
    basealgorithm = CVODE_BDF() 

#=
    Tests on vanilla examples. 
=#
    # Test on vanilla ODE. 
    prob = createsimpleODEproblem(c_tilde, sigma_tilde, mu_tilde, x, T, rho)
    sol = solve(prob, basealgorithm)
    @assert(issorted(sol[end]))

    # Test on vanille PADE. 
    probDAE = createsimpleDAEproblem(c_tilde, sigma_tilde, mu_tilde, x, T, rho)
    solDAE = solve(probDAE, IDA())
    plot(solDAE, vars=1:plotevery:M)
    @show(issorted(solDAE[end][1:M]))

    #Check they are "reasonably" close
    @show norm(sol[1] - solDAE[1][1:M])
    @show norm(sol[end] - solDAE[end][1:M])

#=
    Test for nonuniform. 
=#

    # Create the nonuniform grid. 
    x_add=linspace(x_min,x_max,10)
    x_comb=unique(sort([x;x_add]))
    M_comb=size(x_comb,1)

    # Test it for the ODE. 
    prob_nonuni = createsimpleODEproblem(c_tilde, sigma_tilde, mu_tilde, x_comb, T, rho)
    sol_nonuni = solve(prob_nonuni, basealgorithm)
    plot(sol, vars=1:plotevery:M)
    @assert(issorted(sol_nonuni[end]))

    # Test closeness to uniform. 
    @show norm(sol[1]-sol_nonuni[1]) # check whether close to uniform solution
    @show norm(sol[end]-sol_nonuni[end])