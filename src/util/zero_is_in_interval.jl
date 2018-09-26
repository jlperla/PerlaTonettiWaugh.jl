function zero_is_in_interval(f::Function, lb::Float64, ub::Float64; trial_count = 100)
    zero_trials = range(lb, stop = ub, length = trial_count) # form solution candidates, uniformly gridded
    find_f_zero = zero_trial -> try find_zero(f, zero_trial) catch; Inf end # find solution for f'(t) = 0 (if not Inf)
    zeros = find_f_zero.(zero_trials) # find zeros from the solution candidates
    return any((zeros .> lb) .* (zeros .< ub)) # check if there is a solution in a candidate
end