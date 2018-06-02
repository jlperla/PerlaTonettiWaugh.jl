# To access with v0.6
You will need to get the file in your path.  Download the repository to a `/projects` instead of your normal `Documents` directory or home directory, then
```
push!(LOAD_PATH, joinpath(homedir(), "projects"))
```
This can be added to your julia startup, for example
