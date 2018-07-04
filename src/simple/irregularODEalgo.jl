#Create DiffEq Problem for solving as a system of ODE using nonuniform grid
function createsimplenonunifromODEproblem(c_tilde, sigma_tilde, mu_tilde, x, M::Int64, T::Float64, rho::Float64)
    x, L_1_plus, L_2  = irregulardiffusionoperators(x, M) #Discretize the operator

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