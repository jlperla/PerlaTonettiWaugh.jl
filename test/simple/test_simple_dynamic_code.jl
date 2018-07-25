# Global constants. 
    solve(x,y)=Sundials.solve(x,y)
    # For the test grid. 
    x_min = 0.0
    x_max = 5.0 
    M = 500 
    x = linspace(x_min, x_max, M) # Since we only use the grid now. 
    # For the time boundary. 
    T = 10.0
    # Individual model parameters. 
    π(t, x) = exp(x)*(1.0+0.0*t);
    σ_n = 0.02
    α = 2.1
    ζ_n = 14.5
    r_n = 0.05
    ζ(t, x) = ζ_n + 0.0*x + 0.0*t
    r(t, x) = r_n + 0.0*x + 0.0*t
    σ(t, x) = σ_n + 0.0*x + 0.0*t
    g(t, x) = g_analytic + 0.0*x + 0.0*t;
    γ=0.005;
    # Basic params_pi and settings objects. 
    params_n= @with_kw (γ=0.005, σ=0.02, α=2.1, r=0.05, ζ=14.5)
    settings = @with_kw (z = z,g = g, T = T)
    params_pi = @with_kw (γ=γ, σ=σ, π=π, ζ=ζ, r = r, α=α)
    basealgorithm = CVODE_BDF() 

#= 
    Solutions. 
=#
    # solve for numerical g_T as test
    res=stationary_numerical_simple(params_n(),z)
    g_analytic=res.g
    v_analytic=res.v
    @show v_analytic

    # 1. test starting from g_analytic, v_T are close to analytical result
    prob = create_dynamic_ODE(params_pi(),settings())
    sol = solve(prob, basealgorithm)
    # plot(sol, vars=1:plotevery:M)
    @test issorted(sol[end])
    #@show sol[end]


    # 2. test whether ODE result close to analytic
    tol=1e-8
    @test norm(v_analytic-sol[end])<tol

    # If start at some time varying function of g(t) (change g grid), notice we shouldn't expect this is close to v_analytic

    g_g(t, x) = g_analytic + 0.0*x + 0.01*t;
    settings_g = settings(g = g_g)

prob_g = create_dynamic_ODE(params_pi(), settings_g())
sol_g = solve(prob_g, basealgorithm)
#@show sol_g[end]
@test issorted(sol_g[end])

# 3. Test non uniform grid z

z_comb= unique([linspace(z_min, 1.0, 500)' linspace(1.0, z_max, 201)'])

settings_ir = @with_kw (z = z_comb,g = g, T = T)

prob_ir = create_dynamic_ODE(params_pi(),settings_ir())
sol_ir = solve(prob_ir, basealgorithm)
#@show sol_ir[end]
@test issorted(sol_ir[end])

# interpolate uniform onto non uniform grid
sol_int=LinInterp(z, sol[end])
@show norm(sol_int.(z_comb)-sol_ir[end])
@show norm(sol_ir[end,end]-sol[end,end])

# If comparing with numerical using z_comb grid (g_analytic will be different)
res_comb=stationary_numerical_simple(params_n(), z_comb)
g_acomb=res_comb.g
v_acomb=res_comb.v
g_comb(t, x) = g_acomb + 0.0*x + 0.0*t;

settings_ir2 = @with_kw (z = z_comb,g = g_comb, T = T)

prob_ir2 = create_dynamic_ODE(params_pi(),settings_ir2())
sol_ir2 = solve(prob_ir2, basealgorithm)

@show norm(v_acomb-sol_ir2[end])
# 7. compare sol_ir2 which is chaning both z and g, with the interpolated sol_ir
@show norm(sol_ir2[end]-sol_ir[end])

# 4. check the uniform operator,, use same g_analytic

settings_uni1 = @with_kw (z = z,g = g, T = T)

prob_uni1 = create_dynamic_ODE(params_pi(),settings_uni1())
sol_uni1 = solve(prob_uni1, basealgorithm)

@test norm(sol[end]-sol_uni1[end])<1e-8

# 5. check the uniform operator, adding points
M_uni=701;
z_uni=linspace(z_min,z_max,M_uni)

settings_uni2 = @with_kw (z = z_uni,g = g, T = T)
prob_uni2 = create_dynamic_ODE(params_pi(), settings_uni2())
sol_uni2 = solve(prob_uni2, basealgorithm)

@show norm(sol_uni2[end,end]-sol_uni1[end,end])

# 6. test results from irregular grid with irregular but more dense grid

z_comb3= unique([linspace(z_min, 1.0, 500)' linspace(1.0,2.0, 120)' linspace(2.0, z_max, 81)'])

settings_ir3 = @with_kw (z = z_comb3,g = g, T = T)

prob_ir3 = create_dynamic_ODE(params_pi(),settings_ir3())
sol_ir3 = solve(prob_ir3, basealgorithm)

# interpolate uniform onto non uniform grid
sol_int3=LinInterp(z_comb3, sol_ir3[end])
@show norm(sol_int.(z_comb)-sol_ir[end])
@show norm(sol_ir3[end,end]-sol_ir[end,end]) #at bottom
@show norm(sol_ir3[1,end]-sol_ir[1,end]) # at top

@show g_analytic
@show g_acomb

# 8. check sol_ir2 with sol, which change both z_grid and g_grid


sol_int2=LinInterp(z_comb, sol_ir2[end])
@show norm(sol_int2.(z)-sol[end])
