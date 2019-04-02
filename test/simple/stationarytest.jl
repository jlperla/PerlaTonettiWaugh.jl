# Generate default parameters.
simple_params = @with_kw (μ = 0.0048, υ = 0.02, θ = 2.1, r = 0.05, ζ = 14.5, ξ = 1.0, π = x -> 1)

# Algebraic Tests
results_algebraic = stationary_algebraic_simple(simple_params());
@test results_algebraic.g ≈ 0.020727753490265812; # from updating the eq in stationary
@test results_algebraic.ν ≈ 1.797254207263066;
@test results_algebraic.v(0) ≈ 34.58676260411018;
@test results_algebraic.v(2)*exp(2) ≈ 164.54095232228548; # correct for rescaling
@test results_algebraic.v(5)*exp(5) ≈ 3298.0717492785516;

# Grid Tests
grids = [
    # old grids 
    unique([range(0.0, stop = 1.0, length = 500)' range(1.0, stop = 5.0, length = 201)']),
    unique([range(0.0, stop = 1.0, length = 1000)' range(1.0, stop = 2.0, length = 11)' range(2.0, stop = 5.0, length = 20)']),
    unique([range(0.0, stop = 1.0, length = 300)' range(1.0, stop = 2.0, length = 50)' range(2.0, stop = 7.0, length = 50)']),
    # new grids 
    range(0., stop = 5., length = 100),
    range(0., stop = 5., length = 300),
    range(0., stop = 5., length = 700),
    unique([range(0., stop = 1., length = 400)' range(1., stop = 5., length = 400)'])
]

gs = [
    0.019903560646320612,
    0.020084770046224866,
    0.020440266069327622,
    0.01664042289147539,
    0.01664042289147539,
    0.019533972685068623,
    0.019862262497972627
]

for el in zip(grids, gs),  
    results = stationary_numerical_simple(simple_params(), el[1])
    @show results
    @test results.g ≈ el[2] # regression
    @show results.g - results_algebraic.g # accuracy 
end 
