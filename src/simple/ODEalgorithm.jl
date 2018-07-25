#Create DiffEq Problem for solving as a system of ODE
function createsimpleODEproblem(c_tilde, sigma_tilde, mu_tilde, x::AbstractArray, T::Float64, rho::Float64)
    x, L_1_minus, L_1_plus, L_2  = diffusionoperators(x) #Discretize the operator

    #Check upwind direction
    bothpos = minimum(mu_tilde.(T, x)) >= 0.0 && minimum(mu_tilde.(0.0, x)) >= 0.0
    bothneg = minimum(mu_tilde.(T, x)) <= 0.0 && minimum(mu_tilde.(0.0, x)) <= 0.0
    @assert bothpos || bothneg

    # Dispatch L_1 based on direction initially
    if all(mu_tilde.(0.0, x) .<= 0)
        L_1 = L_1_minus
    elseif all(mu_tilde.(0.0, x) .>= 0)
        L_1 = L_1_plus
    else 
        error("Not weakly positive or negative") # Not strictly required. 
    end

    p = @NT(L_1 = L_1, L_2 = L_2, x = x, rho = rho, mu_tilde = mu_tilde, sigma_tilde = sigma_tilde, c_tilde = c_tilde, T = T) #Named tuple for parameters.

    #Calculating the stationary solution,
    L_T = rho*I - Diagonal(mu_tilde.(T, x)) * L_1 - Diagonal(sigma_tilde.(T, x).^2/2.0) * L_2
    u_T = L_T \ c_tilde.(T, x)

    @assert(issorted(u_T)) #We are only solving versions that are increasing for now

    function f(du,u,p,t)
        L = (p.rho*I  - Diagonal(p.mu_tilde.(t, x)) * p.L_1 - Diagonal(p.sigma_tilde.(t, p.x).^2/2.0) * p.L_2)
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