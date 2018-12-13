@testset "entry_residuals tests" begin 
    @testset "entry_residuals with default Ω" begin
        # Define common objects. 
        params = parameter_defaults()

        settings = settings_defaults()
        settings = merge(settings, (global_transition_penalty_coefficient = 1.0, ))
        z_grid = settings.z
        M = length(z_grid)
    
        d_0 = params.d_0
        d_T = params.d_T
        params_0 = merge(params, (d = d_T,)) # parameters to be used at t = 0
        params_T = merge(params, (d = d_0,)) # parameters to be used at t = T
    
        # solve for stationary solution at t = 0
        stationary_sol_0 = stationary_numerical(params_0, z_grid) # solution at t = 0
        stationary_sol = stationary_numerical(params_T, z_grid) # solution at t = T
    
        Ω_0 = stationary_sol_0.Ω;
        Ω_T = stationary_sol.Ω;
        settings = merge(settings, (params_T = params_T, stationary_sol_T = stationary_sol, Ω_0 = Ω_0))
    
        E_nodes_interior = zeros(length(settings.transition_x0))
        residuals = residuals_given_E_nodes_interior(E_nodes_interior, settings)
        @test mean(residuals) ≈ 0.0 atol = 1e-3
    end

    @testset "entry_residuals with constant Ω" begin
        # Define common objects. 
        params = parameter_defaults()

        settings = settings_defaults()
        settings = merge(settings, (global_transition_penalty_coefficient = 1.0, ))
        z_grid = settings.z
        M = length(z_grid)
    
        d_0 = params.d_0
        d_T = params.d_T
        params_0 = merge(params, (d = d_T,)) # parameters to be used at t = 0
        params_T = merge(params, (d = d_T,)) # parameters to be used at t = T
    
        # solve for stationary solution at t = 0
        stationary_sol_0 = stationary_numerical(params_0, z_grid) # solution at t = 0
        stationary_sol = stationary_numerical(params_T, z_grid) # solution at t = T
    
        Ω_0 = stationary_sol_0.Ω;
        Ω_T = stationary_sol.Ω;
        settings = merge(settings, (params_T = params_T, stationary_sol_T = stationary_sol, Ω_0 = Ω_0))
    
        E_nodes_interior = settings.transition_x0        
        residuals = residuals_given_E_nodes_interior(E_nodes_interior, settings)
        @test mean(residuals) ≈ 0.0 atol = 1e-3
    end
end 