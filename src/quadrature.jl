# New method. 
function ω_weights(z, α, ξ)
     #= 
        Notation: z consists of M points z_1, z_2, ..., z_M.  
        Formula: For the interior points, the weights are given by [f(z_k)/2 * (Δ_k + Δ_{k+1})]
        where Δ_k = z_k - z_{k-1}. For the boundary points, the weights are given by [f(z_1)/2 * Δ_2]
        and [f(z_M)/2 * Δ_M]
    =# 
    Δ = diff(z)
    M = length(z)
    prepend!(Δ, NaN) # To keep the indexing straight. Now, Δ[2] = Δ_2 = z_2 - z_1. And NaN will throw an error if we try to use it.
    # Define the distribution.
    @assert z[1] == 0.0 # Check that our minimum is 0.0 
    z_bar = z[end]
    # Get the vector of probability weights. 
    f_vec = (α * exp.(z * (ξ - α)))/(1 - exp(-α*z_bar))
    # Turn these into weights.
    interiorWeights = [f_vec[i]/2 * (Δ[i] + Δ[i+1]) for i = 2:M-1]
    return [f_vec[1]/2 * Δ[2]; interiorWeights; f_vec[M]/2 * Δ[M]]
end 
