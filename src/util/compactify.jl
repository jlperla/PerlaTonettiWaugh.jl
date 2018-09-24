abstract type TransformerType end
struct UnitCircle <: TransformerType end
struct Logistic <: TransformerType end

# define logit/logistic functions
logistic(x::Float64) = 1.0 / (1.0 + exp(-x))
logit(x::Float64) = log(x / (1.0 - x))

# define Compactifier
struct Compactifier{T<:TransformerType}
    lb::Float64
    ub::Float64
    k::Float64 # half of the length between lb and ub
    a::Float64 # lb + k
    Compactifier(lb::Float64, ub::Float64) = Compactifier(lb::Float64, ub::Float64, UnitCircle())
    Compactifier(lb::Float64, ub::Float64, transformer_type::TransformerType) = lb >= ub ? error("ub should be strictly greater than lb.") : new{typeof(transformer_type)}(lb, ub, (ub - lb)/2, (ub + lb)/2)
end
(f::Compactifier{T})(x) where T <: UnitCircle = f.a + f.k * x / sqrt(1 + x*x)
(f::Compactifier{T})(x) where T <: Logistic = f.lb + (f.ub - f.lb) * logistic(x) 

# define Decompactifier
struct Decompactifier
    compactifier::Compactifier
    Decompactifier(compactifier::Compactifier) = new(compactifier)
end

decompactify_approximately(f::Compactifier{T}, y::Float64) where T <: UnitCircle = (y - f.a) / sqrt(1 - ((y-f.a)/f.k)^2) 
decompactify_approximately(f::Compactifier{T}, y::Float64) where T <: Logistic = logit((y - f.lb) / (f.ub - f.lb)) 
(f::Decompactifier)(y) = decompactify_approximately(f.compactifier, y)