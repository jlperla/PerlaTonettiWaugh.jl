@testset "solve_system" begin
        as = (a -> [a; (a+1)]).(0:10)

        for a in as
            h(x,a) = (x - a) # residual measures how much x is different from a
            x0 = [5.0; 5.0]
            @test solve_system(x -> h(x,a), x0; lb = fill(0.0, length(x0))) ≈ a
        end 
end
