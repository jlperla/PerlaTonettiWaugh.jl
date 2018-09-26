g = x -> x
@test !zero_is_in_interval(g, 0.1, 0.2)
@test !zero_is_in_interval(g, -0.2, -0.1)
@test zero_is_in_interval(g, -0.2, 0.2)
@test zero_is_in_interval(g, -0.5, 0.5)

g = x -> x^2 - 1 # solution is x = 1 and -1
@test !zero_is_in_interval(g, -0.2, 0.2)
@test !zero_is_in_interval(g, -0.5, 0.5)
@test !zero_is_in_interval(g, -2.0, -1.5)
@test !zero_is_in_interval(g, 1.5, 2.0)
@test zero_is_in_interval(g, -1.5, 1.5)
@test zero_is_in_interval(g, -1.5, 0.0)
@test zero_is_in_interval(g, 0.0, 1.5)