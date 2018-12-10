# returns the m-vector zero of a multivariate function h: R^m -> R^n, using nlopt
# constraints_fg! takes (h, x, jacobian_t) and assigns jacobian_t and h[:] given current x.
function find_zero(h, x0; lb = fill(0.0, length(x0)), constraints_fg! = nothing, constraints_tol = fill(1e-8, length(x0)))
    function f(x)
        resids = h(x)
        return sum(resids .* resids)
    end

    function g!(G::Vector, x::Vector)
        ForwardDiff.gradient!(G, f, x)
    end

    function fg!(x::Vector, grad::Vector)
        if length(grad) > 0 # gradient of f(x)
            g!(grad, x)
        end
        f(x)
    end

    # define the optimization problem
    opt = Opt(:LD_SLSQP, length(x0)) # 3 indicates the length of `x`
    lower_bounds!(opt, lb) # find `x` above 0
    min_objective!(opt, fg!) # specifies that optimization problem is on minimization
    if (constraints_fg! != nothing) # add constraints_fg! if needed
        inequality_constraint!(opt, constraints_fg!, constraints_tol)
    end

    xtol_rel!(opt, -Inf)
    xtol_abs!(opt, -Inf)
    ftol_rel!(opt, -Inf)
    ftol_abs!(opt, -Inf)

    # solve the optimization problem
    (minf,minx,ret) = NLopt.optimize(opt, x0)

    return minx
end