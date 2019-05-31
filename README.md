[![Build Status](https://travis-ci.com/jlperla/PerlaTonettiWaugh.jl.svg?token=G6ge79qYLosYiRGJBp1G&branch=master)](https://travis-ci.com/jlperla/PerlaTonettiWaugh.jl)

# Overview

This repository has the complete code for the steady-state and transition dynamics of Perla, Tonetti, and Waugh [Equilibrium Technology Diffusion, Trade, and Growth](http://christophertonetti.com/files/papers/PerlaTonettiWaugh_DiffusionTradeAndGrowth.pdf)

The [derivation document](/docs/numerical_algorithm.pdf) has the complete set of equations implemented for the model, where all equation numbers in the code refer to this document.  General code and derivations for upwind finite difference methods are in the [SimpleDifferentialOperators.jl](https://github.com/QuantEcon/SimpleDifferentialOperators.jl) package with [detailed derivations](https://github.com/QuantEcon/SimpleDifferentialOperators.jl/releases/download/dev/discretized-differential-operator-derivation.pdf).

As in the derivation, the code has a "warmup" model without trade or monopolistic competition to understand transition dynamics with this sort of growth model, and for experimenting with the DAE and finite-difference discretization methods.

## Installation and Use

1. Follow the instructions to [install Julia and Jupyter](https://lectures.quantecon.org/jl/getting_started.html)

2. Open the Julia REPL (see the documentation above) and then install the package (by entering package mode) with

    ```julia
    ] add https://github.com/jlperla/PerlaTonettiWaugh.jl.git
    ```

3. There are several ways you can run the notebooks after installation

    Using the built-in Jupyter is straightforward.  In the Julia terminal
    ```julia
    using PerlaTonettiWaugh, IJulia
    notebook(detached=true, dir=dirname(dirname(pathof(PerlaTonettiWaugh))))
    ```

    Alternatively, to use a separate Jupyter installation you may have installed with Anaconda,
    ```julia
    using PerlaTonettiWaugh
    cd(dirname(dirname(pathof(PerlaTonettiWaugh))))
    ; jupyter lab
    ```
    where the last step runs your `jupyter lab` in the shell.

    **Note** In either case, the first time the `using` it will be very slow

4. Within the main directory, the notebook `full_transition_dynamics.ipynb` has code to solve the model in steady state and on a transition path between steady states.
   - This will cache the output in the `data` folder, with a hash for that particular combination of settings and parameters.
   - The associated `.json` in the data file gives the hashed settings and parameters
5. There is also the `simple_transition_dynamics.ipynb` notebook to solve the simple variation of the model described in the notes.

**NOTE:** When using the notebooks for the first time, it will be very slow as the package and its dependencies are all compiled.
