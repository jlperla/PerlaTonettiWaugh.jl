# perform compactification in reasonable ranges
logistic(x::Float64) = 1.0 / (1.0 + exp(-x))
logit(x::Float64) = log(x / (1.0 - x))

compactify(x::Float64, lb::Float64, ub::Float64) = lb + (ub - lb) * logistic(x) # TODO: compactify
decompactify(x::Float64, lb::Float64, ub::Float64) = logit((x - lb) / (ub - lb)) # TODO: decompactify

get_compactifier(lb::Float64, ub::Float64) = x -> compactify(x, lb, ub)
get_decompactifier(lb::Float64, ub::Float64) = x -> decompactify(x, lb, ub)