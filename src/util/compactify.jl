struct Compactifier
    lb::Float64
    ub::Float64
    k::Float64 # half of the length between lb and ub
    a::Float64 # lb + k
    Compactifier(lb::Float64, ub::Float64) = lb >= ub ? error("ub should be strictly greater than lb.") : new(lb, ub, (ub - lb)/2, (ub + lb)/2)
end
(f::Compactifier)(x) = f.a + f.k * x / sqrt(1 + x*x)
decompactify_approximately(f::Compactifier, y::Float64) = (y - f.a) / sqrt(1 - ((y-f.a)/f.k)^2) # decompactify with some numerical errors