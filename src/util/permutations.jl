# given an array of N vectors, find all possible permutations of N-length vector whose 
# ith element is an element in the ith vector, using recursion
function get_all_permutations(vs, N)
    if (N < 2) return vs end
    rest_ps = get_all_permutations(N != 2 ? vs[2:end] : vs[end], N - 1)
    ps = []
    for i in vs[1]
        for rest_p in rest_ps
            push!(ps, [i; rest_p])
        end
    end
    return ps
end