# New method, accepts a Distribution
function irregulartrapezoidweights(x, d)
    #= 
        Notation: z consists of M points z_1, z_2, ..., z_M.  
        Formula: For the interior points, the weights are given by [f(z_k)/2 * (Δ_k + Δ_{k+1})]
        where Δ_k = z_k - z_{k-1}. For the boundary points, the weights are given by [f(z_1)/2 * Δ_2]
        and [f(z_M)/2 * Δ_M]
    =# 
    
    # Check input and get grid of probability densities.
    ndims(x) == 1 ? true : error("The supplied  state grid is not one-dimensional.")
    f_vec = pdf.(Ref(d), x)
    @assert x[end] ≈ d.upper atol = 1e-5 && x[1] ≈ d.lower atol = 1e-5 # Check that we get most of the probability mass. 

    # Calculate results and return 
    M = length(f_vec)
    Δ = diff(x)
    prepend!(Δ, NaN) # To keep the indexing straight. Now, Δ[2] = Δ_2 = z_2 - z_1. And NaN will throw an error if we try to use it.

    interiorWeights = [f_vec[i]/2 * (Δ[i] + Δ[i+1]) for i = 2:M-1]
    return [f_vec[1]/2 * Δ[2]; interiorWeights; f_vec[M]/2 * Δ[M]]
end 
