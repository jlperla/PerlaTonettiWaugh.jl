module PerlaTonettiWaugh

# package code goes here
export diffusionoperators, createODEproblem, createDAEproblem
using DifferentialEquations, NamedTuples, BenchmarkTools, BandedMatrices, Sundials

#These would be generated automatically by DiffEqOperators.jl
#These are the differential operators for positive drift 0upwind finite differences
function diffusionoperators(x_min, x_max, M)
    x = linspace(x_min, x_max, M)
    Δ = x[2]-x[1]

    dl_1 = zeros(M-1)
    d_1 = -1 * ones(M)
    d_1[end] = 0
    du_1 = ones(M-1)
    L_1_plus = Tridiagonal(dl_1, d_1, du_1)/Δ

    dl_2 = ones(M-1)
    d_2 = -2 * ones(M)
    d_2[1] = -1
    d_2[end] = -1
    du_2 = ones(M-1)
    L_2 = Tridiagonal(dl_2, d_2, du_2)/(Δ^2)

    #BandedMatrix are much faster, probably because of better specializations in the composition
    return (x, BandedMatrix(L_1_plus, (1, 1)), BandedMatrix(L_2, (1, 1))) #The (1,1) are the off-diagonal bandwidths
end

#Create DiffEq Problem for solving
function createODEproblem(c_tilde, sigma_tilde, mu_tilde, x_min::Float64, x_max::Float64, M::Int64, T::Float64, rho::Float64)
    x, L_1_plus, L_2  = diffusionoperators(x_min, x_max, M) #Discretize the operator

    p = @NT(L_1_plus = L_1_plus, L_2 = L_2, x = x, rho = rho, mu_tilde = mu_tilde, sigma_tilde = sigma_tilde, c_tilde = c_tilde, T = T) #Named tuple for parameters.

    #Check upwind direction
    @assert minimum(mu_tilde.(T, x)) >= 0.0
    @assert minimum(mu_tilde.(0.0, x)) >= 0.0

    #Calculating the stationary solution,
    L_T = rho*I - Diagonal(mu_tilde.(T, x)) * L_1_plus - Diagonal(sigma_tilde.(T, x).^2/2.0) * L_2
    u_T = L_T \ c_tilde.(T, x)

    @assert(issorted(u_T)) #We are only solving versions that are increasing for now

    function f(du,u,p,t)
        L = (p.rho*I  - Diagonal(p.mu_tilde.(t, x)) * p.L_1_plus - Diagonal(p.sigma_tilde.(t, p.x).^2/2.0) * p.L_2)
        A_mul_B!(du,L,u)
        du .-= p.c_tilde.(t, p.x)
    end

    #Checks on the residual
    du_T = zeros(u_T)
    f(du_T, u_T, p, T)
    #@show norm(du_T)
    @assert norm(du_T) < 1.0E-10

    return ODEProblem(f, u_T, (T, 0.0), p)
end

function createDAEproblem(c_tilde, sigma_tilde, mu_tilde, x_min, x_max, M, T, rho)
    x, L_1_plus, L_2  = diffusionoperators(x_min, x_max, M) #Discretize the operator

    p = @NT(L_1_plus = L_1_plus, L_2 = L_2, x = x, rho = rho, mu_tilde = mu_tilde, sigma_tilde = sigma_tilde, c_tilde = c_tilde, M = M) #Named tuple for parameters.

    #Check upwind direction
    @assert minimum(mu_tilde.(T, x)) >= 0.0
    @assert minimum(mu_tilde.(0.0, x)) >= 0.0

    #Calculating the stationary solution,
    L_T = rho*I - Diagonal(mu_tilde.(T, x)) * L_1_plus - Diagonal(sigma_tilde.(T, x).^2/2.0) * L_2
    u_T = L_T \ c_tilde.(T, x)

    u_ex_T = [u_T; 1.0] #Closed form for the trivial linear function we are adding

    @assert(issorted(u_T)) #We are only solving versions that are increasing for now


    function f(resid,du_ex,u_ex,p,t)
        L = (p.rho*I  - Diagonal(p.mu_tilde.(t, x)) * p.L_1_plus - Diagonal(p.sigma_tilde.(t, p.x).^2/2.0) * p.L_2)
        u = u_ex[1:p.M]
        resid[1:M] .= L * u_ex[1:p.M] - p.c_tilde.(t, p.x)
        resid[M+1] .= u_ex[M+1] - 1.0
        resid .-= du_ex
    end

    #Should Verifying that the initial condition is solid
    resid_T = zeros(u_ex_T) #preallocation
    du_ex_T = zeros(u_ex_T)
    f(resid_T, du_ex_T, u_ex_T, p, T)
    @show norm(resid_T)
    @assert norm(resid_T) < 1.0E-10

    tspan = (T, 0.0)
    return DAEProblem(f, zeros(u_ex_T), u_ex_T, tspan, differential_vars = [trues(u_T); false], p)
end


end # module
