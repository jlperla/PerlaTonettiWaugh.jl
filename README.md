[![Build Status](https://travis-ci.com/jlperla/PerlaTonettiWaugh.jl.svg?token=G6ge79qYLosYiRGJBp1G&branch=master)](https://travis-ci.com/jlperla/PerlaTonettiWaugh.jl)

# To Install
- On v0.6, this will require your github ID to install.
```julia
Pkg.clone("http://github.com/jlperla/PerlaTonettiWaugh.jl.git")
```

- To see where it is installed,
```julia
julia> Pkg.dir("PerlaTonettiWaugh")
"C:\\Users\\jlperla\\.julia\\v0.6\\PerlaTonettiWaugh"
```
- Then consider dragging and dropping this into Github Desktop or GitKraken in order to make changes.

# To Use
At that point, in the REPL you can
```
using PerlaTonettiWaugh
```

To run the full regression test,
```
Pkg.test("PerlaTonettiWaugh")
```
or run the underlying tests by loading up the file and `ctrl-enter` or `shift-enter` in Atom or VS Code.
