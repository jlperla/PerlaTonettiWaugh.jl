# Generate default parameters.
params = parameters_simple()

# Algebraic Tests
results_algebraic = stationary_algebraic_simple(params, settings_simple());
@test results_algebraic.g ≈ 0.020727753490265812; # from updating the eq in stationary
@test results_algebraic.ν ≈ 1.797254207263066;
@test results_algebraic.v(0) ≈ 34.58676260411018;
@test results_algebraic.v(2)*exp(2) ≈ 164.54095232228548; # correct for rescaling
@test results_algebraic.v(5)*exp(5) ≈ 3298.0717492785516;

# Grid Tests
grids = [
    unique([range(0.0, stop = 1.0, length = 500)' range(1.0, stop = 5.0, length = 201)']),
    unique([range(0.0, stop = 1.0, length = 1000)' range(1.0, stop = 2.0, length = 11)' range(2.0, stop = 5.0, length = 20)']),
    unique([range(0.0, stop = 1.0, length = 300)' range(1.0, stop = 2.0, length = 50)' range(2.0, stop = 7.0, length = 50)']),
    unique([range(0., stop = 0.1, length = 20)' range(0.1, stop = 1., length = 20)' range(1., stop = 5., length = 20)'])
]

gs = [
    0.019903560646302203,
    0.02008477004627781,
    0.020440266030220186,
    0.02080684556396779 # this is the most accurate, and uses the least gridpoints
]

@test var(gs) < 1e-6 # consistency

for i in 1:length(grids)
    results = stationary_numerical_simple(params, settings_simple(z_ex = grids[i]))
    @test results.g ≈ gs[i] # invariance
end