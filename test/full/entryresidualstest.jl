@testset "entry_residuals tests" begin
    @testset "entry_residuals with default Ω" begin
        # Define common objects.
        params = parameter_defaults()
        settings = settings_defaults()

        E_nodes_interior = zeros(length(settings.transition_x0))
        residuals = weighted_residuals(E_nodes_interior, params, settings)
        @test mean(residuals) ≈ 0.0 atol = 1e-2
    end

    @testset "entry_residuals with constant Ω" begin
        # Define common objects.
        params = parameter_defaults()
        params = merge(params, (d_0 = params.d_T,)) # gives us constant Ω
        settings = settings_defaults()
 
        E_nodes_interior = settings.transition_x0
        residuals = weighted_residuals(E_nodes_interior, params, settings)
        @test mean(residuals) ≈ 0.0 atol = 1e-2
    end
end
