#Create DiffEq Problem for solving as a system of DAE using nonuniform grid
function createsimplenonuniformDAEproblem(c_tilde, sigma_tilde, mu_tilde, x, M::Int64, T::Float64, rho::Float64)
    x, L_1_minus, L_1_plus, L_2  = irregulardiffusionoperators(x, M) #Discretize the operator

    #Check upwind direction

    if all(mu_tilde.(0.0, x) .>= 0)
        p = @NT(L_1 = L_1_plus, L_2 = L_2, x = x, rho = rho, mu_tilde = mu_tilde, sigma_tilde = sigma_tilde, c_tilde = c_tilde, T = T, M=M) #Named tuple for parameters.
    elseif all(mu_tilde.(0.0, x) .<= 0)
        p = @NT(L_1 = L_1_minus, L_2 = L_2, x = x, rho = rho, mu_tilde = mu_tilde, sigma_tilde = sigma_tilde, c_tilde = c_tilde, T = T, M=M) #Named tuple for parameters.
    else 
        error("Not weakly positive or negative") # Not strictly necessary, but good to have redundancy here.
    end
    #Calculating the stationary solution,
    L_T = rho*I - Diagonal(mu_tilde.(T, x)) * p.L_1 - Diagonal(sigma_tilde.(T, x).^2/2.0) * L_2
    u_T = L_T \ c_tilde.(T, x)

    u_ex_T = [u_T; 1.0] #Closed form for the trivial linear function we are adding

    @assert(issorted(u_T)) #We are only solving versions that are increasing for now


    function f(resid,du_ex,u_ex,p,t)
        L = (p.rho*I  - Diagonal(p.mu_tilde.(t, x)) * p.L_1 - Diagonal(p.sigma_tilde.(t, p.x).^2/2.0) * p.L_2)
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