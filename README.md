[![Build Status](https://travis-ci.com/jlperla/PerlaTonettiWaugh.jl.svg?token=G6ge79qYLosYiRGJBp1G&branch=master)](https://travis-ci.com/jlperla/PerlaTonettiWaugh.jl)

# Overview

## Installation and Use

1. Follow the instructions to [install Julia](https://lectures.quantecon.org/jl/getting_started.html)

2. Install the package, open the Julia REPL (see the documentation above) and then run

```julia
using Pkg
pkg"dev https://github.com/jlperla/PerlaTonettiWaugh.jl.git"
```

The folder will be installed in `~/.julia/dev/PerlaTonettiWaugh` on mac/linux and in somewhere like `C:\Users\USERNAME\.julia\dev\PerlaTonettiWaugh` on Windows

If you are having trouble finding it, then (in the REPL) run
```julia
joinpath(DEPOT_PATH[1], "dev", "PerlaTonettiWaugh")
```

3. Startup a terminal or Windows Powershell and navigate to that directory.  Then type
```bash
jupyter lab
```
to open a browser with the Jupyter notebooks

4. Within the main directory, the notebook `solve-transition.ipynb` has code to solve the model in steady state and on a transition path between steady states.
   * After opening the notebook, `Ctrl-A` to select the whole notebook and then `Shift-Enter` to run it, or just shift enter between cells
    * The first time it runs, **it will a long time** since it needs to install package dependencies.  Afterwards, it will be faster when you open the notebook.  However, it will still take a few minutes each time you run it the first time (and a few seconds thereafter).
5. There is also the `solve-simple-transition.ipynb` notebook to solve the simple variation of the model described in the notes.

## Weave Instructions

After installing Weave, run:

```
using IJulia, notebook() # Click yes to the install 
using Conda
Conda.add("nbconvert")
```

To generate the notebooks, run 

```
julia generate.jl
```
