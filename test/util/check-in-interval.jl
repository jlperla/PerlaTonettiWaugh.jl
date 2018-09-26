@testset "zero_is_in_interval" begin 
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

    g = x -> (x-1)^2 - 1 # solution is x = 2 and 0
    @test !zero_is_in_interval(g, 0.8, 1.2)
    @test !zero_is_in_interval(g, 0.5, 1.5)
    @test !zero_is_in_interval(g, -1.0, -0.5)
    @test !zero_is_in_interval(g, 2.5, 3.0)
    @test zero_is_in_interval(g, -2.5, 2.5)
    @test zero_is_in_interval(g, -2.5, 1.0)
    @test zero_is_in_interval(g, 1.0, 2.5)

    g = x -> (x+1)^2 - 1 # solution is x = 1 and -1
    @test !zero_is_in_interval(g, -1.2, -0.8)
    @test !zero_is_in_interval(g, -1.5, -0.5)
    @test !zero_is_in_interval(g, -3.0, -2.5)
    @test !zero_is_in_interval(g, 0.5, 1.0)
    @test zero_is_in_interval(g, -0.5, 0.5)
    @test zero_is_in_interval(g, -2.5, -1.0)
    @test zero_is_in_interval(g, -1.0, 0.5)
end

@testset "is_positive_in_interval" begin 
    g = x -> x
    @test is_positive_in_interval(g, 0.1, 0.2)
    @test !is_positive_in_interval(g, -0.2, -0.1)
    @test !is_positive_in_interval(g, -0.2, 0.2)
    @test !is_positive_in_interval(g, -0.5, 0.5)
    
    g = x -> x^2 - 1 # solution is x = 1 and -1
    @test !is_positive_in_interval(g, -0.2, 0.2)
    @test !is_positive_in_interval(g, -0.5, 0.5)
    @test is_positive_in_interval(g, -2.0, -1.5)
    @test is_positive_in_interval(g, 1.5, 2.0)
    @test !is_positive_in_interval(g, -1.5, 1.5)
    @test !is_positive_in_interval(g, -1.5, 0.0)
    @test !is_positive_in_interval(g, 0.0, 1.5)

    
    g = x -> (x-1)^2 - 1 # solution is x = 2 and 0
    @test !is_positive_in_interval(g, 0.8, 1.2)
    @test !is_positive_in_interval(g, 0.5, 1.5)
    @test is_positive_in_interval(g, -1.0, -0.5)
    @test is_positive_in_interval(g, 2.5, 3.0)
    @test !is_positive_in_interval(g, -2.5, 2.5)
    @test !is_positive_in_interval(g, -2.5, 1.0)
    @test !is_positive_in_interval(g, 1.0, 2.5)

    g = x -> (x+1)^2 - 1 # solution is x = 1 and -1
    @test !is_positive_in_interval(g, -1.2, -0.8)
    @test !is_positive_in_interval(g, -1.5, -0.5)
    @test is_positive_in_interval(g, -3.0, -2.5)
    @test is_positive_in_interval(g, 0.5, 1.0)
    @test !is_positive_in_interval(g, -0.5, 0.5)
    @test !is_positive_in_interval(g, -2.5, -1.0)
    @test !is_positive_in_interval(g, -1.0, 0.5)
end