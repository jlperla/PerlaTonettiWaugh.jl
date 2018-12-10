@testset "find_zero" begin
    @testset "find_zero without nonlinear constraints" begin
        as = (a -> [a; (a+1)]).(0:10)

        for a in as
            h(x,a) = (x - a) # residual measures how much x is different from a
            x0 = [5.0; 5.0]
            @test find_zero(x -> h(x,a), x0; lb = fill(0.0, length(x0))) ≈ a
        end
    end

    @testset "find_zero with nonlinear constraints" begin
        function constraints_fg!(h, x, jacobian_t)
            M = length(x)
            # A is a matrix whose ith row has 1 in ith col and -1 in (i+1)th col
            # so ith element of A*x imposes x[i] <= x[i+1]
            # note that this imposes x[end] <= 0 as well as the last row is [0; 0; ...; 0; 1].
            A = LinearAlgebra.Tridiagonal(zeros(M-1),ones(M),-ones(M-1))
            if length(jacobian_t) > 0 # transpose of the Jacobian matrix
                jacobian_t[:] = A'
            end
            h[:] = A*x
        end
        
        as = (a -> [-a; -a+1; 0.0]).(2:10)
        
        for a in as
            h(x,a) = (x - a)
            x0 = [-8.0; -7.0; 0.0]
            @test find_zero(x -> h(x,a), x0; lb = fill(-20.0, length(x0)), constraints_fg! = constraints_fg!) ≈ a 
        end
    end
end