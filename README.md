[![Build Status](https://travis-ci.com/jlperla/PerlaTonettiWaugh.jl.svg?token=G6ge79qYLosYiRGJBp1G&branch=master)](https://travis-ci.com/jlperla/PerlaTonettiWaugh.jl)

# Overview

Can put some text here using the link to [the paper](http://jesseperla.com/papers/perla_tonetti_waugh.pdf), 

## Installation and Use

1. Follow the instructions to [install Julia](https://lectures.quantecon.org/jl/getting_started.html#Installing-Julia-and-Dependencies)

2. Install the package by running 

```julia 
using Pkg 
pkg"dev https://github.com/jlperla/PerlaTonettiWaugh.jl.git" 
```

Consider dragging and dropping the folder into Github Desktop or GitKraken in order to make changes. 

The folder will be installed in `~/.julia/dev/PerlaTonettiWaugh`. 

The .julia folder is somewhere like `C:\Users\USERNAME\.julia\` on Windows, and `~/.julia` on mac/Linux.

3. Open the above folder in a Jupyter notebook (i.e., go there in a terminal and run `jupyter lab`). If you can't find it, run 

```julia
joinpath(DEPOT_PATH[1], "dev", "PerlaTonettiWaugh")
```

to get the path. 

4. The notebook `solve-transition.ipynb` has code to solve the model in steady state and on a transition path between steady states. 

5. At the end, there is a call to the `solve_continuation` method, which returns a relatively robust solution by smoothly changing the `d_0` parameter from a known steady-state.
