# returns the m-vector zero of a multivariate function h: R^m -> R^n, using nlopt
# constraints_fg! takes (h, x, jacobian_t) and assigns jacobian_t and h[:] given current x.
function find_zero(h, x0; lb = fill(0.0, length(x0)), ub = fill(10e8, length(x0)), 
                    constraints_fg! = nothing, constraints_tol = fill(1e-8, length(x0)),
                    autodiff = :forward)
    function f(x)
        resids = h(x)
        return sum(resids .* resids)
    end

    # define the optimization problem
    opt = (constraints_fg! == nothing) ? Opt(:LD_LBFGS, length(x0)) : Opt(:LN_COBYLA, length(x0)) # 3 indicates the length of `x`
    f_opt = NLoptAdapter(f, x0, autodiff)
    min_objective!(opt, f_opt) # specifies that optimization problem is on minimization

    if (lb != nothing)    
        lower_bounds!(opt, lb) # find `x` above lb
    end
    if (ub != nothing)    
        upper_bounds!(opt, ub) # find `x` below ub
    end

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

struct NLoptAdapter{T} <: Function where T <: AbstractObjective
    nlsolver_base::T
end

# implement fg!; note that the order is reversed
(adapter::NLoptAdapter)(x, df) = adapter.nlsolver_base.fdf(df, x)
(adapter::NLoptAdapter)(result, x, jacobian_transpose) = adapter.nlsolver_base.fdf(result, jacobian_transpose', x)

# constructors
NLoptAdapter(f, x, autodiff = :forward) = NLoptAdapter(OnceDifferentiable(f, x, autodiff = autodiff))
NLoptAdapter(f!, x::Vector, F::Vector, autodiff = :forward) = NLoptAdapter(OnceDifferentiable(f!, x, F, autodiff = autodiff))