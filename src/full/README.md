### README for code associated with solving the transition dynamics and steady state of of [Equilibrium Technology Diffusion, Trade, and Growth](https://christophertonetti.com/files/papers/PerlaTonettiWaugh_DiffusionTradeAndGrowth.pdf) by Perla, Tonetti, and Waugh (AER 2020) 

#### Julia Source

The files in this package are installed and executed through installation in the top-level [README.md](../../README.md). The files are:

  - [``params.jl``](params.jl) organizes parameters, settings, and initial conditions.
    - All of these may be swapped out using the named tuples.
    - To speed up loading pre-solved versions of the model, we store a cache in the `/data` folder where the name is defined by `model_cachename(...)` function in this file.  The cachename is calculated by hashing all parameters and settings.
    - `load_parameters(..)` takes a CSV file with calibrated parameters (generated from the `/src/calibration` etc. ) and creates the necessary structure for the Julia files.
- [``static.jl``](static.jl) directly maps equations from the paper for the static equilibrium calculations from the parameters and intermediate values.
- [``stationary.jl``](stationary.jl) calculates the stationary equilibrium and associated welfare analysis.
  - In particular `stationary_algebraic(parameters, settings)`  solves the model as the system of equations specified in the main paper, and using the `static.jl` functions
  - `stationary_numerical(parameters, settings)` solves for the equilibrium using upwind-finite difference methods for the ODEs rather than solving as a system of equations.  This is primarily used as the initial condition for the transition dynamics (which require ODEs).
  - `steady_state_from_g` is used by the `total_derivative` the welfare analysis
- [``dynamic.jl``](dynamic.jl) calculates the transition dynamics of the equilibrium between the two steady-states
  - The main entry-point is `solve_transition(parameters, settings)` which calculates the two steady states for the system of DAEs.
  - That, in turn iterates on the stock of varieties, `Ω` and solves the system of DAEs conditional on that sequence with `solve_dynamics(...)`.  Convergence occurs when the entry residual is minimized, while the adoption decision is encapsulated in the DAE.
  - Finally, `prepare_results` generates a dataframe from the results.