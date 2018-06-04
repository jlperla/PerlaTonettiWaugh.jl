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
- Then consider dragging and dropping this into Github Desktop in order to make changes.

# To Use
At that point, you can go
```
using PerlaTonettiWaugh
```

To run the full regression test,
```
Pkg.test("PerlaTonettiWaugh")
```
